import math
import random
from prettytable import PrettyTable
import matplotlib.pyplot as plt


class Param:
    h, dis, w, d, J = None, None, None, None, None
    alpha = []

    def __init__(self, h1, dis1, w1, d1, J1, alpha1):
        self.alpha.clear()
        self.h = h1
        self.dis = dis1
        self.w = w1
        self.d = d1
        self.J = J1
        self.alpha = alpha1


    def __lt__(self, other):
        return self.J < other.J


class Signal:
    result_list = []

    signal_vec = []
    noise_vec = []
    x_vec = []
    filtered_vec = []
    alpha_vec = []

    params_vec = []

    points_vec = []

    x_min = 0
    x_max = math.pi
    K = 100
    P = 0.95
    L = 10
    e = 0.01
    r = 3
    M = int((r - 1) / 2)
    N = int(math.log(1 - P) / math.log(1 - (e / (x_max - x_min))))

    def __init__(self):
        self.signal_vec.clear()
        self.noise_vec.clear()
        self.filtered_vec.clear()
        self.result_list.clear()
        self.x_vec.clear()
        self.alpha_vec.clear()
        self.params_vec.clear()
        self.points_vec.clear()

        for i in range(self.K):
            self.x_vec.append(self.x_min + i * (self.x_max - self.x_min) / self.K)
            self.signal_vec.append(math.sin(self.x_vec[i]) + 0.5)

        for i in range(self.K):
            self.noise_vec.append(self.signal_vec[i] + random.uniform(-0.25, 0.25))

        self.result_list.append(self.signal_vec)
        self.result_list.append(self.noise_vec)

    def generate_alpha(self):
        temp_alpha = []
        for i in range(0, self.r):
            temp_alpha.append(0)
        temp_alpha[self.M] = float(format(random.uniform(0, 1), '.4f'))

        sum = 0.0
        for m in range(2, int(self.M + 1)):
            for s in range(m, self.r - m):
                sum += temp_alpha[s]
            temp_alpha[m - 1] = float(format(0.5 * random.uniform(0, 1 - sum), '.4f'))
            temp_alpha[self.r - m] = float(format(0.5 * random.uniform(0, 1 - sum), '4f'))

        sum = 0.0
        for s in range(1, self.r - 1):
            sum += temp_alpha[s]
        temp_alpha[0] = float(format(0.5 * (1 - sum), '.4f'))
        temp_alpha[self.r - 1] = float(format(0.5 * (1 - sum), '.4f'))

        sum = 0.0
        for i in temp_alpha:
            sum += i
        for i in temp_alpha:
            i /= sum

        return temp_alpha

    def make_filter(self, k, temp_alpha):
        mult = 1
        for j in range(k - self.M, k + self.M):
            mult *= pow(self.noise_vec[j], temp_alpha[j + self.M - k])
        return mult

    def make_filtered_signal(self):
        self.alpha_vec = self.generate_alpha()
        for k in range(self.M, self.K - self.M):
            self.filtered_vec.append(self.make_filter(k, self.alpha_vec))
        self.result_list.append(self.filtered_vec)

    def random_search_method(self, el):
        self.r = el
        self.M = int((self.r - 1) / 2)

        for i in range(0, self.L + 1):
            temp_lambda = i / self.L

            min_params = []
            for j in range(0, self.N):
                temp_alpha = self.generate_alpha()

                temp_w = 0.0
                for k in range(1, self.K - 2 * self.M):
                    temp_w += abs(self.make_filter(k, temp_alpha) - self.make_filter(k - 1, temp_alpha))

                temp_d = 0.0
                for k in range(0, self.K - 2 * self.M):
                    temp_d += abs(self.make_filter(k, temp_alpha) - self.noise_vec[k])
                temp_d /= self.K

                temp_J = temp_lambda * temp_w + (1 - temp_lambda) * temp_d
                distance = abs(temp_w) + abs(temp_d)
                min_params.append(Param(temp_lambda, distance, temp_w, temp_d, temp_J, temp_alpha))

            self.params_vec.append(min(min_params))

        self.make_filtered_signal()

        x = PrettyTable()
        x.field_names = ["h", "dis", "alpha", "w", "d"]
        for el in self.params_vec:
            x.add_row([el.h, "%.4f" % el.dis, el.alpha, "%.4f" % el.w, "%.4f" % el.d])
        print(x)
        x.clear()

        min_el = self.params_vec[0].dis
        min__ = self.params_vec[0]
        for el in self.params_vec:
            if el.dis < min_el:
                min__ = el
                min_el = el.dis

        y = PrettyTable()
        y.field_names = ["h*", "J", "w", "d"]
        y.add_row([min__.h, "%.4f" % min__.J, "%.4f" % min__.w, "%.4f" % min__.d])
        print(y)
        y.clear()

        for el in self.params_vec:
            self.points_vec.append([el.w, el.d])

        # построеноие графиков
        plt.plot(self.result_list[0], label='f(x) = sin(x) + 0.5')
        plt.plot(self.result_list[1], label='noise')
        plt.plot(self.result_list[2], label='filtering')
        plt.legend(title=None, loc='lower center')
        plt.title('Functions')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.grid(True)
        plt.show()

        for point in self.points_vec:
            plt.plot(point[0], point[1], 'o')
        plt.show()

if __name__ == '__main__':
    signal_v1 = Signal()
    signal_v1.random_search_method(3)

    signal_v2 = Signal()
    signal_v2.random_search_method(5)
