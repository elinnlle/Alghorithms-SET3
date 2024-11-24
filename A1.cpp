#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

// Функция для проверки, попадает ли точка внутрь круга
bool isInCircle(double cx, double cy, double r, double x, double y) {
    return (x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r;
}

// Функция Монте-Карло для оценки площади пересечения трёх кругов
void areaMonteCarlo(const vector<double>& a, const vector<double>& b, const vector<double>& c, bool narrowBounds) {
    double x_min, x_max, y_min, y_max;

    // Определение границ области
    if (narrowBounds) {
        x_min = max({a[0] - a[2], b[0] - b[2], c[0] - c[2]});
        x_max = min({a[0] + a[2], b[0] + b[2], c[0] + c[2]});
        y_min = max({a[1] - a[2], b[1] - b[2], c[1] - c[2]});
        y_max = min({a[1] + a[2], b[1] + b[2], c[1] + c[2]});
    } else {
        x_min = min({a[0] - a[2], b[0] - b[2], c[0] - c[2]});
        x_max = max({a[0] + a[2], b[0] + b[2], c[0] + c[2]});
        y_min = min({a[1] - a[2], b[1] - b[2], c[1] - c[2]});
        y_max = max({a[1] + a[2], b[1] + b[2], c[1] + c[2]});
    }

    // Генератор случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist_x(x_min, x_max);
    uniform_real_distribution<> dist_y(y_min, y_max);

    // Истинная площадь (если известна аналитически)
    double exactArea = 0.25 * M_PI + 1.25 * asin(0.8) - 1.0;

    // Заголовок таблицы
    cout << setw(10) << "N"
         << setw(27) << "Площадь"
         << setw(31) << "Погрешность" << endl;

    // Увеличиваем число точек в последовательных экспериментах
    for (int numPoints = 100; numPoints <= 100000; numPoints += 500) {
        int pointsInside = 0;

        for (int i = 0; i < numPoints; ++i) {
            double x = dist_x(gen);
            double y = dist_y(gen);

            // Проверяем, находится ли точка внутри всех трёх кругов
            if (isInCircle(a[0], a[1], a[2], x, y) &&
                isInCircle(b[0], b[1], b[2], x, y) &&
                isInCircle(c[0], c[1], c[2], x, y)) {
                pointsInside++;
            }
        }

        // Вычисляем оценку площади для данного числа точек
        double areaRect = (x_max - x_min) * (y_max - y_min);
        double estimatedArea = (static_cast<double>(pointsInside) / numPoints) * areaRect;

        // Вычисляем погрешность (если известна точная площадь)
        double error = fabs(estimatedArea - exactArea) / exactArea;

        // Выводим результаты
        cout << setw(10) << numPoints
             << setw(20) << fixed << setprecision(6) << estimatedArea
             << setw(20) << fixed << setprecision(6) << error << endl;
    }
}

int main() {
    // Параметры трёх окружностей (x, y, r)
    vector<double> a = {1.0, 1.0, 1.0};
    vector<double> b = {1.5, 2.0, sqrt(5) / 2};
    vector<double> c = {2.0, 1.5, sqrt(5) / 2};

    // Запускаем алгоритм Монте-Карло
    cout << "\nРезультаты для узкой области:\n";
    areaMonteCarlo(a, b, c, true);

    cout << "\nРезультаты для широкой области:\n";
    areaMonteCarlo(a, b, c, false);

    return 0;
}
