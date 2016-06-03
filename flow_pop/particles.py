import numpy as np

class Simulation(object):

    def __init__(self, Lx=1., Ly=1., interaction_length=0.02, num_particles=300,
                 num_populations = 2):
        self.Lx = Lx
        self.Ly = Ly
        self.interaction_length = interaction_length

        self.num_particles = num_particles
        self.num_populations = num_populations

        self.Nx = np.int(self.Lx / self.interaction_length)
        self.Ny = np.int(self.Ly / self.interaction_length)

        starting_num = self.num_particles / self.num_populations

        # Positions of each type of particle
        self.positions = np.random.rand(starting_num, 2, self.num_populations)

        scale = np.array([self.Lx, self.Ly])

        self.positions *= scale

        self.grid_positions = np.int32(self.positions / interaction_length)

        self.grid = np.zeros((self.Nx, self.Ny, self.num_populations), dtype=np.int32)

        for cur_pop in range(self.num_populations):
            xy = self.grid_positions[:, :, cur_pop]
            for cur_xy in xy:
                self.grid[cur_xy[0], cur_xy[1]] += 1