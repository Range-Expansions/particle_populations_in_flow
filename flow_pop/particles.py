import numpy as np

class Simulation_2d(object):

    def __init__(self, Lx=1., Ly=1., interaction_length=0.02, num_particles=300,
                 num_populations = 2, dt = 0.1):

        self.Lx = Lx
        self.Ly = Ly
        self.L = np.array([Lx, Ly])
        self.dt = dt

        self.interaction_length = interaction_length

        self.num_particles = num_particles
        self.num_populations = num_populations

        self.num_bins = self.L/interaction_length

        starting_num = num_particles / num_populations

    class Particle(object):
        def __init__(self, simulation, pop_type, position, grid_point, D=1.0):

            self.simulation = simulation

            self.pop_type = np.int32(pop_type)
            self.position = np.float32(position)
            self.grid_point = np.int32(grid_point)

            self.D = D

        def move(self):

            grid[self.grid_point[0], self.grid_point[1], self.pop_type] -= 1

            self.position += np.sqrt(2 * self.D * dt) * np.random.randn(2)

            # Deal with moving out of the system...bounce back
            Lx = L[0]
            Ly = L[1]

            x = self.position[0]
            y = self.position[1]

            dx = 0

            if (x < 0):
                x = np.abs(x)
            if (x > Lx):
                x = Lx - (x - Lx)

            if (y < 0):
                y = np.abs(y)
            if (y > Ly):
                y = Ly - (y - Ly)

            self.position[0] = x
            self.position[1] = y

            self.grid_point = np.int32(self.position / interaction_length)

            grid[self.grid_point[0], self.grid_point[1], self.pop_type] += 1

        def react(self):