import numpy as np
import weakref


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

        self.num_bins = np.int32(self.L/interaction_length)

        # Create particles randomly for now
        starting_num = num_particles / num_populations

        self.particle_list = []
        for cur_pop in range(num_populations):
            for i in range(starting_num):
                cur_position = np.random.rand(2) * self.L

                cur_grid = np.int32(cur_position / interaction_length)

                new_particle = Particle(self, cur_pop, cur_position, cur_grid)

                self.particle_list.append(new_particle)

        # Setup the grid
        self.grid = np.zeros((self.num_bins[0], self.num_bins[1], self.num_populations), dtype=np.int32)
        for cur_particle in self.particle_list:
            xy = cur_particle.grid_point
            pop_num = cur_particle.pop_type

            self.grid[xy[0], xy[1], pop_num] += 1

    def react(self):
        """Right now, simple concentration-based growth"""

        particles_to_keep = []

        for cur_particle in self.particle_list:
            x = cur_particle.grid_point[0]
            y = cur_particle.grid_point[1]

            if cur_particle.pop_type != self.num_populations: # Last type is the concentration field

                num_c = self.grid[x, y, self.num_populations]

                prob = num_c * cur_particle.k * self.dt
                rand = np.random.rand()

                if rand < prob: # React!
                    particles_to_keep.append(cur_particle)

                    new_particle = cur_particle.birth()
                    # Keep the old particle
                    particles_to_keep.append(new_particle)

            elif cur_particle.pop_type == self.num_populations:

                for i in range(self.num_populations - 1): # Loop over all possible 2-pair interactions
                    cur_n = self.grid[x, y, i]
                    prob = cur_n * cur_particle.k * self.dt
                    rand = np.random.rand()

                    if rand < prob: # Die
                        pass
                    else: # Keep the particle around
                        particles_to_keep.append(cur_particle)

        # Appropriately add and delete particles to the simulation after the loop
        self.particle_list = particles_to_keep

    def move(self):
        for cur_particle in self.particle_list:
            cur_particle.move()

    def run(self, num_iterations):
        for i in range(num_iterations):
            self.move()
            self.react()

class Particle(object):
    def __init__(self, simulation, pop_type, position, grid_point, D=1.0, k = 1.0):

        self.sim = weakref.proxy(simulation)

        self.pop_type = np.int32(pop_type)
        self.position = np.float32(position)
        self.grid_point = np.int32(grid_point)

        self.D = D
        self.k = k

    def move(self):

        sim = self.sim

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] -= 1

        self.position += np.sqrt(2 * self.D * sim.dt) * np.random.randn(2)

        # Deal with moving out of the system...bounce back
        Lx = sim.L[0]
        Ly = sim.L[1]

        x = self.position[0]
        y = self.position[1]

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

        self.grid_point = np.int32(self.position / sim.interaction_length)

        sim.grid[self.grid_point[0], self.grid_point[1], self.pop_type] += 1

    def birth(self):
        return Particle(self.sim, self.pop_type, self.position, self.grid_point,
                        D=self.D, k=self.k)