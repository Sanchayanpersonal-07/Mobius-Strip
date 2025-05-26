# mobius_strip_model.py

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D # Not strictly needed for recent matplotlib versions

class MobiusStrip:
    """
    Models a Mobius strip using parametric equations and computes geometric properties.
    """
    def __init__(self, R, w, n):
        """
        Initializes the MobiusStrip.

        Args:
            R (float): Radius of the center loop of the strip.
            w (float): Width of the strip.
            n (int): Resolution parameter. This primarily defines the number of points
                     along the length (u-direction) of the strip. The number of points
                     across the width (v-direction) will be derived from this.
        """
        if R <= 0:
            raise ValueError("Radius R must be positive.")
        if w <= 0:
            raise ValueError("Width w must be positive.")
        if w >= 2 * R:
            print("Warning: Width w might be too large relative to R, potentially leading to self-intersections not handled by this simple model.")
        if n <= 1:
            raise ValueError("Resolution n must be greater than 1.")

        self.R = R
        self.w = w
        self.resolution_param = n

        # Define n_u (points along length) and n_v (points across width)
        self.n_u = n
        self.n_v = max(3, int(self.n_u / 10.0)) # Ensure at least 3 points for width

        # Placeholder for computed properties
        self.X = None
        self.Y = None
        self.Z = None
        self.U_grid = None
        self.V_grid = None

        self._surface_area = None
        self._edge_length = None

        self._generate_mesh()

    def _parametric_eq(self, u, v):
        """
        Parametric equations for the Mobius strip.
        u in [0, 2*pi], v in [-w/2, w/2]
        """
        x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        z = v * np.sin(u / 2)
        return x, y, z

    def _generate_mesh(self):
        """
        Generates the 3D mesh points (X, Y, Z) for the Mobius strip surface.
        The points are stored in 2D arrays self.X, self.Y, self.Z.
        Shape of X, Y, Z will be (n_v, n_u) due to meshgrid's 'xy' indexing.
        """
        u_domain = np.linspace(0, 2 * np.pi, self.n_u, endpoint=True)
        v_domain = np.linspace(-self.w / 2, self.w / 2, self.n_v, endpoint=True)

        self.U_grid, self.V_grid = np.meshgrid(u_domain, v_domain, indexing='xy')
        self.X, self.Y, self.Z = self._parametric_eq(self.U_grid, self.V_grid)

    @property
    def mesh_points_stacked(self):
        """
        Returns the mesh points as a single array of shape (N_points, 3).
        N_points = n_u * n_v.
        """
        return np.stack((self.X.ravel(), self.Y.ravel(), self.Z.ravel()), axis=-1)

    def compute_surface_area(self):
        """
        Numerically approximates the surface area of the Mobius strip
        by summing the areas of small triangular patches on the surface.
        """
        if self._surface_area is not None:
            return self._surface_area

        area = 0.0
        for j in range(self.n_v - 1):
            for i in range(self.n_u - 1):
                p00 = np.array([self.X[j, i],   self.Y[j, i],   self.Z[j, i]])
                p10 = np.array([self.X[j, i+1], self.Y[j, i+1], self.Z[j, i+1]])
                p01 = np.array([self.X[j+1, i], self.Y[j+1, i], self.Z[j+1, i]])
                p11 = np.array([self.X[j+1,i+1],self.Y[j+1,i+1],self.Z[j+1,i+1]])

                vec1_tri1 = p10 - p00
                vec2_tri1 = p01 - p00
                area += 0.5 * np.linalg.norm(np.cross(vec1_tri1, vec2_tri1))

                vec1_tri2 = p11 - p10
                vec2_tri2 = p01 - p10
                area += 0.5 * np.linalg.norm(np.cross(vec1_tri2, vec2_tri2))
        
        self._surface_area = area
        return self._surface_area

    def compute_edge_length(self):
        """
        Numerically approximates the length of the single continuous edge of the Mobius strip.
        """
        if self._edge_length is not None:
            return self._edge_length

        length = 0.0

        edge_pts_top_x = self.X[-1, :] 
        edge_pts_top_y = self.Y[-1, :]
        edge_pts_top_z = self.Z[-1, :]

        edge_pts_bot_x = self.X[0, :]
        edge_pts_bot_y = self.Y[0, :]
        edge_pts_bot_z = self.Z[0, :]

        for i in range(self.n_u - 1):
            p1 = np.array([edge_pts_top_x[i], edge_pts_top_y[i], edge_pts_top_z[i]])
            p2 = np.array([edge_pts_top_x[i+1], edge_pts_top_y[i+1], edge_pts_top_z[i+1]])
            length += np.linalg.norm(p2 - p1)

        p_top_last = np.array([edge_pts_top_x[-1], edge_pts_top_y[-1], edge_pts_top_z[-1]])
        p_bot_first = np.array([edge_pts_bot_x[0], edge_pts_bot_y[0], edge_pts_bot_z[0]])
        length += np.linalg.norm(p_bot_first - p_top_last)

        for i in range(self.n_u - 1):
            p1 = np.array([edge_pts_bot_x[i], edge_pts_bot_y[i], edge_pts_bot_z[i]])
            p2 = np.array([edge_pts_bot_x[i+1], edge_pts_bot_y[i+1], edge_pts_bot_z[i+1]])
            length += np.linalg.norm(p2 - p1)
        
        p_bot_last = np.array([edge_pts_bot_x[-1], edge_pts_bot_y[-1], edge_pts_bot_z[-1]])
        p_top_first = np.array([edge_pts_top_x[0], edge_pts_top_y[0], edge_pts_top_z[0]])
        length += np.linalg.norm(p_top_first - p_bot_last)

        self._edge_length = length
        return self._edge_length

    def plot(self, ax=None, show_edges=True, save_path=None):
        """
        Plots the Mobius strip using Matplotlib.

        Args:
            ax (matplotlib.axes.Axes, optional): Existing 3D axes to plot on.
            show_edges (bool, optional): Whether to draw the edge of the strip.
            save_path (str, optional): If provided, saves the plot to this file path.
        """
        if self.X is None:
            print("Mesh not generated. Call _generate_mesh() or ensure initialization is complete.")
            return

        if ax is None:
            fig = plt.figure(figsize=(12, 9))
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.get_figure()

        ax.plot_surface(self.X, self.Y, self.Z, cmap='viridis', 
                        edgecolor='k', linewidth=0.2, alpha=0.7, 
                        rstride=1, cstride=1)

        if show_edges:
            plot_edge_x = np.concatenate((self.X[-1, :], self.X[0, :], [self.X[-1,0]]))
            plot_edge_y = np.concatenate((self.Y[-1, :], self.Y[0, :], [self.Y[-1,0]]))
            plot_edge_z = np.concatenate((self.Z[-1, :], self.Z[0, :], [self.Z[-1,0]]))
            ax.plot(plot_edge_x, plot_edge_y, plot_edge_z, color='red', linewidth=3, label='Edge Path')

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title(f'Möbius Strip (R={self.R}, w={self.w}, n_u={self.n_u}, n_v={self.n_v})')
        
        all_X = self.X.flatten()
        all_Y = self.Y.flatten()
        all_Z = self.Z.flatten()
        max_range = np.array([all_X.max()-all_X.min(), all_Y.max()-all_Y.min(), all_Z.max()-all_Z.min()]).max() / 2.0
        mid_x = (all_X.max()+all_X.min()) * 0.5
        mid_y = (all_Y.max()+all_Y.min()) * 0.5
        mid_z = (all_Z.max()+all_Z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        if show_edges:
            ax.legend()

        if save_path:
            plt.savefig(save_path)
            print(f"Plot saved to {save_path}")
        
if __name__ == "__main__":
    # Define parameters for the Mobius strip
    R_center = 2.0
    strip_width = 1.0
    resolution = 100

    try:
        mobius = MobiusStrip(R=R_center, w=strip_width, n=resolution)
    except ValueError as e:
        print(f"Error creating Mobius strip: {e}")
        exit()

    print(f"Möbius Strip Parameters:")
    print(f"  Center Radius (R): {mobius.R}")
    print(f"  Strip Width (w): {mobius.w}")
    print(f"  Input Resolution (n): {mobius.resolution_param}")
    print(f"  Mesh Resolution (n_u x n_v): {mobius.n_u} x {mobius.n_v} points")
    print("-" * 30)

    surface_area = mobius.compute_surface_area()
    print(f"Approximate Surface Area: {surface_area:.4f}")

    edge_length = mobius.compute_edge_length()
    print(f"Approximate Edge Length: {edge_length:.4f}")
    print("-" * 30)

    print("Generating plot...")
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    mobius.plot(ax=ax, show_edges=True, save_path="mobius_strip_visualization.png")
    
    plt.tight_layout()
    plt.show()