// МЕТОД ИТЕРАЦИЙ
#include <iostream>
#include <cmath>
//class Tianwen_1 {
//private:
//    double apocenter = 12000000;
//    double periapsis = 265000;
//public:
//    semi_major_axis(double apocenter; double periapsis) {
//
//    }
//}

  long double eccentricity = 0.9567787607;
  const long double PI = 3.14159265358979323846;
  int period = 28080; //c
//long double calculateEccentricAnomaly(long double eccentricity,long double mean_anomaly,long double precision) {
//    long double E0 = mean_anomaly;
//    long double diff = 0.001; // Diff - это переменная, которая представляет разницу между двумя последовательными значениями  в процессе итерации
//    int max_iterations = 10000;
//    int iteration = 0;
//    while (diff > precision && iteration < max_iterations) {
//        long double E1 = mean_anomaly + eccentricity * sin(E0);
//        diff = abs(E1 - E0);
//        E0 = E1;
//        iteration++;
//    }
//    return E0;
//}
//
//int main() {
//   
//    long double eccentricity = 0.9567787607;
//    long double precision = 1e-10; // требуемая точность
//    for (int time = 0; time <= 28080; time++) {
//        long double mean_anomaly = 2 * PI * time / 28080; // средняя аномалия для каждого момента времени
//        long double eccentric_anomaly = calculateEccentricAnomaly(eccentricity, mean_anomaly, precision);
//        std::cout << "Time: " << time << ", Eccentric Anomaly: " << eccentric_anomaly << std::endl;
//    }
//    return 0;
//}

// ПОЛОВИННОЕ ДЕЛЕНИЕ



// функция для вычисления средней аномалии по эксцентрической аномалии
//double meanAnomaly(double eccentricity, double eccentricAnomaly) {
//    return eccentricAnomaly - eccentricity * sin(eccentricAnomaly);
//}
//
//// функция для вычисления эксцентрической аномалии методом половинного деления
//double calculateEccentricAnomaly(double meanAnomaly, double eccentricity, double precision) {
//    // начальные значения интервала [a, b]
//    double a = -2.0 * PI;
//    double b = 3.0 * PI;
//
//    // итерации с уменьшением точности
//    while (fabs(b - a) >= precision) {
//        double c = (a + b) / 2.0;
//        double functionValue = c - eccentricity * sin(c) - meanAnomaly;
//
//        if (functionValue == 0.0)
//            break;
//        else if (functionValue * (a - b) < 0)
//            b = c;
//        else
//            a = c;
//    }
//
//    return (a + b) / 2.0;
//}
//
//int main() {
//   
//
//    // итерации для каждого момента времени от 0 до 28080
//    for (int time = 0; time <= 28080; ++time) {
//        // приведение момента времени к радианам
//        double meanAnomalyRad = (2.0 * PI * time) / period;
//
//        // вычисление эксцентрической аномалии с точностью от 10^(-3) до 10^(-10)
//        for (int i = 3; i <= 10; ++i) {
//            double precision = pow(10, -i);
//            double eccentricAnomaly = calculateEccentricAnomaly(meanAnomalyRad, eccentricity, precision);
//            double calculatedMeanAnomaly = meanAnomaly(eccentricity, eccentricAnomaly);
//
//            std::cout << "time: " << time <<  ", eccentric anomaly:" << eccentricAnomaly
//                << ", calculated mean anomaly: " << calculatedMeanAnomaly << std::endl;
//        }
//    }
//
//    return 0;
//}
// 
// 
// 





  //МЕТОД ЗОЛОТОГО СЕЧЕНИЯ
// Функция для вычисления значения функции F(E) = M + e*sin(E) - E
//double calculateFunction(double M, double e, double E) {
    //return M + e * sin(E) - E;
}

// Функция для вычисления эксцентрической аномалии методом золотого сечения
//double calculateEccentricAnomaly(double M, double e, double accuracy) {
    //const double goldenRatio = (1 + sqrt(5)) / 2; // Золотое сечение

    //double a = 0;
    //double b = 2 * M_PI; // Пределы для эксцентрической аномалии
    //double c = b - (b - a) / goldenRatio;
   // double d = a + (b - a) / goldenRatio;

    //while (abs(b - a) > 2 * accuracy) {
        //if (calculateFunction(M, e, c) * calculateFunction(M, e, d) < 0) {
           // b = d;
       // }
        //else {
            //a = c;
       // }

       // c = b - (b - a) / goldenRatio;
       // d = a + (b - a) / goldenRatio;
    //}

    // Найденная эксцентрическая аномалия
    //double E = (b + a) / 2;
    //return E;
//}

//int main() {
    
    //double accuracy = 1e-10; // Точность вычислений

    //for (double time = 0; time <= period; time += 1) {
        //double M = 2 * M_PI * time / 28080; // Средняя аномалия
        //double E = calculateEccentricAnomaly(M, eccentricity, accuracy);

        //std::cout << "Time: " << time << ", Eccentric Anomaly: " << E << std::endl;
    //}

    //return 0;
//}
 



  // МЕТОД НЬЮТОНА

 
