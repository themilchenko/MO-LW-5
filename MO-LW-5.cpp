/* Вариант задания (14):                      */
/*      Среднее геометрическое, Манхэттенская */

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#define PI 3.14159265

/*структтурирование данных из сводной таблицы*/
struct Params
{
    double h;
    double dis;
    std::vector<double> alpha_vec;
    double w;
    double d;
    double J;

    Params() = default;
    Params(double h1, double dis1, const std::vector<double> alpha, double w1, double d1, double j) :
        h(h1), dis(dis1), w(w1), d(d1), J(j)
    {
        for (const double& i : alpha)
            alpha_vec.push_back(i);
    }

    Params(const Params& other) : h(other.h), dis(other.dis), w(other.w), d(other.d), J(other.J)
    {
        for (const double& i : other.alpha_vec)
            alpha_vec.push_back(i);
    }

    Params& operator=(const Params& other)
    {
        if (this != &other)
        {
            if (!alpha_vec.empty())
                alpha_vec.clear();
            h = other.h;
            dis = other.dis;
            w = other.w;
            d = other.d;
            J = other.J;
            alpha_vec.resize(other.alpha_vec.size());
            for (int i = 0; i < other.alpha_vec.size(); ++i)
                alpha_vec[i] = other.alpha_vec[i];
        }
        return *this;
    }
};

class Signals
{
    std::vector<std::pair<double, double>> signal_vec;              /*вектор исходного сигнала*/
    std::vector<std::pair<double, double>> noise_vec;               /*вектор шума*/

    std::vector<double> filtered_vec;                               /*отфильтрованный шум*/

    std::vector<Params> params_vec;                                 /*данные из таблицы*/

    /*===============*/
    /*~~~константы~~~*/
    /*===============*/
    const double x_min_ = 0;
    const double x_max_ = PI;
    const int K_ = 100;

    const double P_ = 0.95;
    const int L_ = 10;
    const double epsilon = 0.01;

    double r_ = 3;
    double M_ = (r_ - 1) / 2;
    const int N_ = static_cast<int>((std::log(1 - P_)) / (std::log(1 - (epsilon / (x_max_ - x_min_)))));

public:
    /*инициализация сигнала и его шума*/
    Signals()
    {
        for (int i = 0; i < K_; ++i)
        {
            std::pair<double, double> element;
            element.first = x_min_ + i * (x_max_ - x_min_) / K_;
            element.second = sin(element.first) + 0.5;
            signal_vec.push_back(element);
        }

        std::random_device r;
        std::default_random_engine generator(r());
        std::uniform_real_distribution<double> distribution(-0.25, 0.25);

        for (int i = 0; i < K_; ++i)
        {
            std::pair<double, double> element;
            element.first = signal_vec[i].first;
            double random_element = distribution(generator);
            element.second = signal_vec[i].second + random_element;
            noise_vec.push_back(element);
        }
    }

    /*функция гегерации альфа ветора*/
    std::vector<double> generate_alpha()
    {
        /*создание вектора и инициализация его размерности*/
        std::vector<double> alpha_vec(r_);

        /*определение генерации псевдослучайных чисел*/
        std::random_device ran;
        std::default_random_engine generator(ran());
        std::uniform_real_distribution<double> distribution(0, 100);

        double random_element = distribution(generator) / 100;

        /*по указанному алгоритму задаю значения альфы*/
        alpha_vec[M_] = random_element;

        double sum = 0.0;
        for (size_t m = 2; m < M_ + 1; ++m)
        {
            for (size_t s = m; s <= r_ - m; ++s)
                sum += alpha_vec[s];

            std::uniform_real_distribution<double> new_distribution(0, 1 - sum);
            //generator.seed(std::time(nullptr));
            random_element = new_distribution(generator);

            alpha_vec[m - 1] = 0.5 * random_element;
            alpha_vec[r_ - m] = 0.5 * random_element;
        }

        sum = 0.0;
        for (size_t s = 1; s < r_; ++s)
            sum += alpha_vec[s];
        alpha_vec[0] = 0.5 * (1 - sum);
        alpha_vec[r_ - 1] = 0.5 * (1 - sum);

        /*нормирование вектора*/
        sum = 0.0;
        for (double i : alpha_vec)
            sum += i;
        for (double i : alpha_vec)
            i /= sum;

        return alpha_vec;
    }

    /*печать линии*/
    void print_line_v1(size_t num)
    {
        std::cout << '+';
        for (size_t i = 0; i < 5; ++i)
            std::cout << '-';
        std::cout << '+';
        for (size_t i = 0; i < 8; ++i)
            std::cout << '-';
        std::cout << '+';
        for (size_t i = 0; i < num; ++i)
            std::cout << '-';
        std::cout << '+';
        for (size_t i = 0; i < 8; ++i)
            std::cout << '-';
        std::cout << '+';
        for (size_t i = 0; i < 8; ++i)
            std::cout << '-';
        std::cout << "+\n";
    }

    void print_line_v2()
    {
        std::cout << "+-----+--------+--------+--------+\n";
    }

    /*печать шапки таблицы*/
    void print_title_v1()
    {
        std::cout << "|  h  |  dis   |          alpha           |   w    |   d    |\n";
    }

    void print_title_v2()
    {
        std::cout << "|  h* |   J    |   w    |   d    |\n";
    }

    /*печать данных таблицы*/
    void print_data_v1(double h, double dist, const std::vector<double>& alpha, double w, double d)
    {
        std::cout << "| " << std::setprecision(1) << std::fixed << h << " |" <<
            " " << std::setprecision(4) << std::fixed << dist << " | [";
        for (size_t i = 0; i < alpha.size(); ++i)
        {
            if (i == alpha.size() - 1)
                std::cout << std::setprecision(4) << std::fixed << alpha[i];
            else
                std::cout << std::setprecision(4) << std::fixed << alpha[i] << ", ";
        }
        std::cout << "] | " << std::setprecision(4) << std::fixed << w <<
            " | " << d << " |\n";
    }

    void print_data_v2(const Params& data)
    {
        std::cout << "| " << std::setprecision(1) << std::fixed << data.h <<
            " | " << std::setprecision(4) << std::fixed << data.J <<
            " | " << std::setprecision(4) << std::fixed << data.w <<
            " | " << std::setprecision(4) << std::fixed << data.d << " |\n";
    }

    /*по ф-ле взвешенного скользящего среднего нахождение значения функции очистки сигнала в точке*/
    double make_filter(size_t k, const std::vector<double>& temp_alpha)
    {
        double mult = 1;
        for (int j = k - M_; j < k + M_; ++j)
        {
            mult *= std::powf(noise_vec[j].second, temp_alpha[j + M_ + 1 - k]);
        }
        return mult;
    }

    void make_filtered_signal(const std::vector<double>& alpha_vec)
    {
        for (int k = 0; k < K_; ++k)
            filtered_vec.push_back(make_filter(k, alpha_vec));
    }

    /*метод случайного поиска*/
    void random_search_method(int r)
    {
        r_ = r;
        M_ = (r_ - 1) / 2;

        print_line_v1(26);
        print_title_v1();
        print_line_v1(26);

        for (int i = 0; i < L_ + 1; ++i)
        {
            std::vector<double> best_alpha; /*наиболее подходящий альфа вектор*/

            double temp_lambda = static_cast<double>(i) / L_;
            double min_J = 999999999;  /*задаие максимального значения переменной J для поиска минимума*/
            double min_w = 0;
            double min_d = 0;
            double distance = 0;

            std::vector<std::vector<double>> alpha_mtrx;
            std::vector<Params> temp_params;

            /*цикл с количеством экспериметов N*/
            for (int j = 0; j < N_; ++j)
            {
                std::vector<double> temp_alpha = generate_alpha(); /*генерция */
                alpha_mtrx.push_back(temp_alpha);
                std::vector<double> temp_filtered_vec;

                /*релизация очищенного сигнала*/
                //for (int k = M_ + 1; k < K_ - M_; ++k)
                //    temp_filtered_vec.push_back(make_filter(k, temp_alpha));

                double temp_w = 0.0;
                for (int k = 1; k < K_; ++k)
                    temp_w += fabs(make_filter(k, temp_alpha) - make_filter(k - 1, temp_alpha));

                double temp_d = 0.0;
                for (int k = 0; k < K_; ++k)
                    temp_d += fabs(temp_filtered_vec[k] - noise_vec[k].second);
                temp_d /= K_;

                double temp_J = temp_lambda * temp_w + (1 - temp_lambda) * temp_d;

                if (temp_J < min_J)
                {
                    min_J = temp_J;
                    min_w = temp_w;
                    min_d = temp_d;
                    best_alpha.resize(temp_alpha.size());
                    for (int i = 0; i < temp_alpha.size(); ++i)
                        best_alpha[i] = temp_alpha[i];
                    distance = fabs(min_w) + fabs(min_d);
                }
                
            }
            print_data_v1(temp_lambda, distance, best_alpha, min_w, min_d);                       /*печать полученных минимальных данных*/
            params_vec.push_back(Params(temp_lambda, distance, best_alpha, min_w, min_d, min_J)); /*копирование данных таблицы в отдкльный вектор*/
        }
        print_line_v1(26);

        /*очистка шумового сигнала*/
        make_filtered_signal(data_processing().alpha_vec);
        print_filtered_vec();

        std::cout << std::endl;

        print_line_v2();
        print_title_v2();
        print_line_v2();
        print_data_v2(data_processing());
        print_line_v2();
    }

    Params data_processing()
    {
        Params result_data;
        result_data.dis = 999999999;
        for (int i = 0; i < params_vec.size(); ++i)
            if (params_vec[i].dis < result_data.dis)
                result_data = params_vec[i];
        return result_data;
    }

    /*вывод сигнала и шума для построения графиков*/
    void print()
    {
        std::cout << "SIGNAL VECTOR:\n";
        for (const auto& el : signal_vec)
            std::cout << el.first << ' ' << el.second << '\n';

        std::cout << "NOISE VECTOR:\n";
        for (const auto& el : noise_vec)
            std::cout << el.second << '\n';
    }

    void print_filtered_vec()
    {
        std::cout << "FILTERED VECTOR:\n";
        for (double el : filtered_vec)
            std::cout << el << '\n';
    }

};

int main()
{
    //Signals signal;
    //signal.print();
    //signal.random_search_method(3);

    Signals signal_v2;
    signal_v2.print();
    signal_v2.random_search_method(5);
}


/*
1,0546
0,8875
0,9706
0,9386
1,0889
1,0970
1,1574
1,1823
1,0013
0,9814
1,1812
1,2419
1,1141
1,1562
1,2023
1,2174
1,2983
1,1820
1,1865
1,2532
1,1551
1,3344
1,1473
1,2965
1,1897
1,2525
1,2005
1,2860
1,3062
1,3583
1,4011
1,2708
1,3915
1,3987
1,4256
1,2715
1,3171
1,3376
1,4587
1,4880
1,4123
1,3097
1,4967
1,4237
1,4226
1,4601
1,4060
1,3695
1,4874
1,4485
1,3265
1,3722
1,5087
1,3489
1,3972
1,3209
1,3470
1,3568
1,4021
1,3182
1,4521
1,4056
1,4111
1,3657
1,4572
1,4727
1,4593
1,3848
1,4104
1,4267
1,2720
1,2614
1,3461
1,2224
1,2130
1,3795
1,3910
1,3323
1,3298
1,1360
1,2843
1,2517
1,2487
1,3053
1,2416
1,2447
1,0268
1,1946
1,1636
1,0064
1,2166
0,9776
0,9999
1,0005
1,1135
1,0264
1,0148
1,0142
1,0600
0,7984
*/