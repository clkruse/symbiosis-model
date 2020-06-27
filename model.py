import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

default_params = {
    "sym_load": 1,
    "ros_tolerance_mean": 0.8,
    "ros_tolerance_variance": 0.1,
    "starvation_time": 14
}

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
        self.ros_tolerance_mean = 0.8
        self.ros_tolerance_variance = 0.1
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
        ros_output = self.symbiont.calculate_ros_production(temperature)
        if ros_output > self.ros_tolerance:
            self.sym_load = self.ros_tolerance / ros_output
        self.health += 2 * self.sym_load / self.starvation_time - 1.9 / self.starvation_time

        return self.health



class Symbiont(object):
    def __init__(self):
        self.ros_output = 0
        self.temp_sensitivity = np.random.uniform(0,2)
        self.sensitivity_center = np.random.normal(28, 1)

    def calculate_ros_production(self, temperature):
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
        if coral.health < 0.9 and coral.health > 0:
            new_symbiont = Symbiont()
            coral.symbiont = new_symbiont

        if coral.health <= 0:
            coral.health = 0
            coral.is_living = False








test_sym = Symbiont()
test_coral = Coral(test_sym)
test_sym.sensitivity_center
test_coral.health = 1

history_health = []
history_temp = []
history_sym = []
day = np.random.randint(100000)
for i in range(1000):
    day += 1
    temperature = generate_temperature(29, day)
    history_temp.append(temperature)
    step(test_coral, temperature)
    history_health.append(test_coral.health)
    history_sym.append(test_coral.sym_load)

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


corals = []
num_corals = 10000
for i in range(num_corals):
    sym = Symbiont()
    corals.append(Coral(sym))

history_health = []
history_temp = []
day = np.random.randint(100, 100000)

for i in trange(1000):
    day += 1
    temperature = generate_temperature(29, day)
    history_temp.append(temperature)

    for coral in corals:
        step(coral, temperature)

print("Survival Rate: {0:.1f}%".format(100 * np.sum([coral.health for coral in corals])/len(corals)))

plt.hist([coral.ros_tolerance for coral in corals], 35)
plt.title("Coral Population ROS Tolerance")
plt.show()

plt.hist([coral.symbiont.calculate_ros_production(30) for coral in corals], 35)
plt.title("Symbiont Population ROS Production at 30°C")
plt.show()



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
