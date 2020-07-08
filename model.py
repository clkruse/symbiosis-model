# %% Imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

# %% Define Functions
def generate_temperature(mean_temperature, day):
    day = day % 365
    year_cycle = 4 * np.cos((np.random.uniform(-20, 20) + day) * 2 * np.pi / 365)
    season_cycle = 2 * np.cos(4 * (np.random.uniform(-20, 20) + day) * 2 * np.pi / 365)
    month_cycle = 1 * np.cos(12 * (np.random.uniform(-5, 5) + day) * 2 * np.pi / 365)
    daily_variation = np.random.normal(0, 0.5)
    #daily_variation = 0

    temperature = mean_temperature + (year_cycle + season_cycle + month_cycle) / 3 + daily_variation

    return temperature

class Coral(object):

    def __init__(self, symbiont):
        self.sym_load = 1
        self.ros_tolerance_mean = 5
        self.ros_tolerance_variance = 0.5
        self.ros_tolerance = np.random.normal(self.ros_tolerance_mean, self.ros_tolerance_variance)
        # Define the time it takes a coral to die of starvation
        # Here, I define it as 14 days
        self.starvation_time = 14
        self.time_without_symbiont = 0
        self.health = 1
        self.symbiont = symbiont
        self.is_living = True

    def symbiont_load(self, ros):
        if ros < self.ros_tolerance:
            self.sym_load = 1
        if ros >= self.ros_tolerance:
            self.sym_load = 0
            self.time_without_symbiont += 1


    def compute_coral_health(self, temperature):
        """
        Update a coral's health based on the temperature.
        Computing health also rebalances the the symbiont load to lower the
        ROS output below the coral's ROS threshold.
        """
        ros_production = self.symbiont.calculate_ros_production(temperature)
        ros_output = self.sym_load * ros_production

        if ros_output > self.ros_tolerance:
            self.sym_load = self.ros_tolerance / ros_output
        self.health += 1 / self.starvation_time * (2 * self.sym_load - 1.5)
        return self.health


class Symbiont(object):
    def __init__(self):
        self.temp_sensitivity = np.random.uniform(1.25,2)
        self.sensitivity_center = np.random.normal(28, 1)

    def calculate_ros_production(self, temperature):
        # Revisiting the temperature response curve
        """
        Calculate symbiont ros production at a given temperature
        Function is an exponential based on a randomly assigned mean and sensitivity.
        https://www.desmos.com/calculator/3gpeckfxrj
        """
        ros_output = self.temp_sensitivity ** (temperature - self.sensitivity_center)

        return ros_output

    def calculate_ros_production_sigmoidal(self, temperature):
        # I'm currently making the symbiont ROS production response to temperature a sigmoid as seen here. https://www.desmos.com/calculator/h7yegs3nsa
        # I don't think this is physiologically accurate
        ros_output = 1/2 * (np.tanh(self.temp_sensitivity * \
                                (temperature - self.sensitivity_center)) + 1)
        return ros_output


def step(coral, temperature):
    if coral.is_living == True:
        coral.health = coral.compute_coral_health(temperature)
        if coral.health > 1:
            coral.health = 1
        if coral.symbiont.calculate_ros_production(temperature) > coral.ros_tolerance and coral.health > 0:
            new_symbiont = Symbiont()
            coral.symbiont = new_symbiont

        if coral.health <= 0:
            coral.health = 0
            coral.is_living = False

# %% Test the process step by step

test_sym = Symbiont()
test_coral = Coral(test_sym)
test_coral.health = 1
test_coral.sym_load = 1
print("Symbiont sensitivity center: {:.2f}".format(test_sym.sensitivity_center))
print("Symbiont temperature sensitivity: {:.2f}".format(test_sym.temp_sensitivity))

print("Coral ROS tolerance: {:.2f}".format(test_coral.ros_tolerance))

#print("ROS Production at 28°C: {:.2f}".format(test_coral.symbiont.calculate_ros_production(28)))
#print("ROS Production at 29°C: {:.2f}".format(test_coral.symbiont.calculate_ros_production(29)))
#print("ROS Production at 30°C: {:.2f}".format(test_coral.symbiont.calculate_ros_production(30)))

temperature = 33
current_coral_health = test_coral.health
print("ROS production at {0}°C: {1:.2f}".format(temperature, test_coral.symbiont.calculate_ros_production(temperature)))

print("Change in coral health at {0}°C: {1:.2f}".format(temperature, test_coral.compute_coral_health(temperature) - current_coral_health))

print("Symbiont load at  {0}°C: {1:.2f}".format(temperature, test_coral.sym_load))

delta_healths = []
symbiont_loads = []
temps = []

for temp in range(1000):
    test_coral.health = 1
    test_coral.sym_load = 1
    temp = 28 + temp/200
    temps.append(temp)
    delta_healths.append(test_coral.compute_coral_health(temp) - 1)
    symbiont_loads.append(test_coral.sym_load)
plt.plot(temps, [health * 10 for health in delta_healths])
plt.title("∆ Health and Sym Load vs. temp")
plt.plot(temps, symbiont_loads, c='r')
plt.legend(['∆ Health * 10', 'Sym Load'])
plt.show()


# %% Individual simulation
test_sym = Symbiont()
test_coral = Coral(test_sym)
test_sym.sensitivity_center
test_coral.health = 1

history_health = []
history_temp = []
history_sym = []
history_ros = []
day = np.random.randint(100000)
for i in range(200):
    day += 1
    temperature = generate_temperature(29, day)
    history_temp.append(temperature)
    step(test_coral, temperature)
    history_health.append(test_coral.health)
    history_sym.append(test_coral.sym_load)
    history_ros.append(test_coral.symbiont.sensitivity_center)

plt.plot(history_health)
plt.title('Health')
plt.ylim([0, 1.1])
plt.show()
plt.plot(history_temp, c='r')
plt.title('Temperature')
plt.show()
plt.plot(history_sym, c='r')
plt.title('Symbiont Load')
plt.show()
plt.plot(history_ros, c='g')
plt.title('ROS Mean Sensitivity')
plt.show()



# %% Population Simulation

def print_population_stats(population):
    print("-"*75)

    mean_coral_ros_tolerance = np.mean([coral.ros_tolerance for coral in population])
    print("Coral ROS tolerance: {:.2f}".format(mean_coral_ros_tolerance))

    print("Coral health: {:.2f}".format(np.mean([coral.health for coral in population])))

    mean_ros_sensitivity_center = np.mean([coral.symbiont.sensitivity_center for coral in population])
    print("Symbiont ROS sensitivity center: {:.2f}".format(mean_ros_sensitivity_center))


    mean_ros_sensitivity = np.mean([coral.symbiont.temp_sensitivity for coral in population])
    print("Symbiont ROS sensitivity: {:.2f}".format(mean_ros_sensitivity))

    print("Symbiont ROS production at 30°C: {:.2f}".format(np.mean([coral.symbiont.calculate_ros_production(30) for coral in population])))

    critical_temperature = mean_ros_sensitivity_center - (np.log(1 / mean_coral_ros_tolerance)) / np.log(mean_ros_sensitivity)

    print("Critical temperature: {:.2f}".format(critical_temperature))

    print("-"*75)


corals = []
num_corals = 5000
for i in range(num_corals):
    sym = Symbiont()
    corals.append(Coral(sym))

history_health = []
history_temp = []
day = np.random.randint(100, 100000)

print("Starting Population Statistics")
print_population_stats(corals)


for i in trange(10000):
    day += 1
    temperature = generate_temperature(29, day)
    history_temp.append(temperature)

    for coral in corals:
        if coral.is_living == True:
            step(coral, temperature)

living_corals = [coral for coral in corals if coral.is_living == True]
dead_corals = [coral for coral in corals if coral.is_living == False]

print("*** Survival Rate: {0:.1f}% ***".format(100 * len(living_corals)/(len(living_corals) + len(dead_corals))))

print("Living Population Statistics")
print_population_stats(living_corals)

print("Dead Population Statistics")
print_population_stats(dead_corals)


# %%
"""
plt.hist([coral.ros_tolerance for coral in living_corals], 35)
plt.title("Coral Population ROS Tolerance")
plt.show()

plt.hist([coral.symbiont.calculate_ros_production(30) for coral in living_corals], 35)
plt.title("Symbiont Population ROS Production at 30°C")
plt.show()

plt.hist([coral.symbiont.temp_sensitivity for coral in living_corals], 35)
plt.title("Temperature Sensitivity {0:.2f}".format(np.mean([coral.symbiont.temp_sensitivity for coral in living_corals])))
plt.show()
"""


# %% Temperature plot
temps = []
random_start = np.random.randint(100, 100000)
for day in range(random_start, random_start + 365 * 2):
    temps.append(generate_temperature(29, day))
plt.figure(figsize=(8,6))
plt.plot(temps, linewidth=1)
plt.title("Temperature - Mean of 29")
plt.xlabel('Day')
plt.ylabel('Temperature')
#plt.savefig("Example Temperature Plot - 2 Years.png", bbox_inches="tight", dpi=200)
plt.show()
