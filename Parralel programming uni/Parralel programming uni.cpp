#include <iostream>
#include <random>
using namespace std;

int main() {
    int N;
    double min = 0, max = 100;

    cout << "Select number of equations: ";
    cin >> N;

    // Generating random double in range
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_real_distribution<> distr(min, max); // define the range

    double** matrix = new double*[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];
        for (int j = 0; j < N; j++) {
            matrix[i][j] = distr(gen);
        }
    }
    


    return 0;
}

/*
void Gaussian_elimination(double(&matrix)[2][3]);
void Zero_line_case(double(&matrix)[2][3], const int& count_1, const int& count_2);
double f(const double& a, const double& x, const double& b, const double& y, const double& c);
double g(const double& a, const double& x, const double& b, const double& y, const double& c);
void result(const double(&original_matrix)[2][3], const double& x, const double& y);

int main() {

    // Решение системы из двух линейных уравнений

    // Ввод и обработка данных
    double matrix[2][3];
    cout << "Solving a system of two linear equations\n"
        << "Enter coefficients a, b and c for the first equation:\n";
    cout << "a = ";
    cin >> matrix[0][0];
    cout << "b = ";
    cin >> matrix[0][1];
    cout << "c = ";
    cin >> matrix[0][2];

    cout << "Enter coefficients a, b and c for the second equation:\n";
    cout << "a = ";
    cin >> matrix[1][0];
    cout << "b = ";
    cin >> matrix[1][1];
    cout << "c = ";
    cin >> matrix[1][2];

    cout << "\nIn total:\n";
    cout << "f(x, y) = " << matrix[0][0] << "x + " << matrix[0][1] << "y = "
        << matrix[0][2] << "\n";
    cout << "g(x, y) = " << matrix[1][0] << "x + " << matrix[1][1] << "y = "
        << matrix[1][2] << "\n\n";


    // Находим выражения для x и y с помощью метода Гаусса
    Gaussian_elimination(matrix);
    return 0;
}


// Метод Гаусса
void Gaussian_elimination(double(&matrix)[2][3]) {

    double x, y;

    // Копируем исходные коэф-ты системы, т.к. значения matrix будут изм-ся
    double original_matrix[2][3];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            original_matrix[i][j] = matrix[i][j];
        }
    }


    // Проверяем на наличие нулевых строк
    int count_1 = 0;    // Число нулевых коэф-ов в 1-ой строке
    int count_2 = 0;    // Число нулевых коэф-ов во 2-ой строке

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            if (matrix[i][j] == 0) {
                if (i == 0) { count_1++; }
                else { count_2++; }
            }
        }
    }

    if (count_1 == 3 || count_2 == 3) {
        Zero_line_case(matrix, count_1, count_2);
        return;
    }
    else {

        // Отсекаем случаи отсутствия решений (т.е. если a == 0, b == 0)
        if (matrix[0][0] == 0 && matrix[0][1] == 0) {
            cout << "System has no solution\n";
            return;
        }
        else if (matrix[1][0] == 0 && matrix[1][1] == 0) {
            cout << "System has no solution\n";
            return;
        }


        // Если matrix[0][0] == matrix[1][0] == 0,
        // то проверяем систему на совместность
        if (matrix[0][0] == 0 && matrix[1][0] == 0) {

            // Складываем модифицированную первую строку 
            // со второй, чтобы получить ноль вместо matrix[1][1]
            // matrix[0][0] и matrix[0][1] не могут одновременно быть == 0,
            // т.к. matrix[0][0] == 0, то проверка matrix[0][1] == 0 не нужна
            double modif = -matrix[1][1] / matrix[0][1];
            matrix[1][1] += modif * matrix[0][1];
            matrix[1][2] += modif * matrix[0][2];


            if (matrix[1][1] == 0 && matrix[1][2] == 0) {
                // Если matrix[1][1] == matrix[1][2] == 0 после сложения 
                // строк, то f(x, y) = k*g(x, y), где k - целое
                // Иными словами, вторую строку можно заменить нулями,
                // а это случай count_2 == 3 и count_1 < 3
                Zero_line_case(matrix, 0, 3);
                return;
            }
            else {
                // При невыполнении предыдущего условия 
                // у системы не может быть решений
                cout << "System has no solution\n";
                return;
            }
        }

        // Если matrix[0][0] == 0, меняем строки местами
        if (matrix[0][0] == 0) {
            double switcher;
            for (int j = 0; j < 3; j++) {
                switcher = matrix[0][j];
                matrix[0][j] = matrix[1][j];
                matrix[1][j] = switcher;
            }

            // Обратный ход метода Гаусса
            // matrix[1][1] и matrix[1][0] не могут одновременно быть == 0,
            // т.к. matrix[1][0] == 0, то проверка matrix[1][1] == 0 не нужна
            y = matrix[1][2] / matrix[1][1];
            x = matrix[0][2] - matrix[0][1] * y;

            // Выводим решение и проверяем подстановкой
            result(original_matrix, x, y);
            return;
        }
        else {
            // Складываем модифицированную первую строку 
            // со второй, чтобы получить ноль вместо matrix[1][0]
            double modif = -matrix[1][0] / matrix[0][0];
            matrix[1][0] = 0;   // Сложение здесь не имеет смысла, сразу ноль
            matrix[1][1] += modif * matrix[0][1];
            matrix[1][2] += modif * matrix[0][2];

            // Обратный ход метода Гаусса
            if (matrix[1][1] != 0) {
                y = matrix[1][2] / matrix[1][1];
                x = (matrix[0][2] - matrix[0][1] * y) / matrix[0][0];

                // Выводим решение и проверяем подстановкой
                result(original_matrix, x, y);
                return;
            }
            else if (matrix[1][1] == 0) {
                // Случай 0*x + 0*y = c, c != 0
                cout << "System has no solution\n";
                return;
            }
            else if (matrix[1][1] == 0 && matrix[1][2] == 0) {
                // Если matrix[1][1] == matrix[1][2] == 0 после сложения 
                // строк, то f(x, y) = k*g(x, y), где k - целое
                // Иными словами, вторую строку можно заменить нулями,
                // а это случай count_2 == 3 и count_1 < 3
                Zero_line_case(matrix, 0, 3);
                return;
            }
        }
    }

}

// Случай нулев(ой|ых) строк
void Zero_line_case(double(&matrix)[2][3], const int& count_1, const int& count_2) {

    // Случай нулевых строк
    if (count_1 + count_2 == 6) {
        cout << "x - and y - any numbers\n";
        return;
    }

    // Случай одной нулевой строки
    // Если первая строка нулевая, то меняем строки местами
    if (count_1 == 3) {
        for (int j = 0; j < 3; j++) {
            matrix[0][j] = matrix[1][j];
            matrix[1][j] = 0;
        }
    }

    // Перебираем различные варианты вида ненулевой строки
    if (count_1 == 3 || count_2 == 3) {
        // a != 0, b != 0, c - any number
        if (matrix[0][0] != 0 && matrix[0][1] != 0) {
            cout << "y = " << matrix[0][2] / matrix[0][1]
                << " - " << matrix[0][0] / matrix[0][1] << "x\n"
                << "x - any number\n";
            return;
        }
        // a == 0, b != 0, c - any number
        else if (matrix[0][0] == 0 && matrix[0][1] != 0) {
            cout << "y = " << matrix[0][2] / matrix[0][1] << '\n'
                << "x - any number\n";
            return;
        }
        // a != 0, b == 0, c - any number
        else if (matrix[0][0] != 0 && matrix[0][1] == 0) {
            cout << "x = " << matrix[0][2] / matrix[0][0] << '\n'
                << "y - any number\n";
            return;
        }
        // a == 0, b == 0, c != 0
        else if (matrix[0][0] == 0 && matrix[0][1] == 0) {
            cout << "System has no solution\n";
            return;
        }
    }
}

// Функции системы
double f(const double& a, const double& x, const double& b, const double& y, const double& c) {
    return a * x + b * y - c;
}

double g(const double& a, const double& x, const double& b, const double& y, const double& c) {
    return a * x + b * y - c;
}

// Вывод результата работы метода на экран
void result(const double(&original_matrix)[2][3], const double& x, const double& y) {
    cout << "x = " << x << ", y = " << y << '\n'

        << "f(x, y) = " << f(original_matrix[0][0], x,
            original_matrix[0][1], y,
            original_matrix[0][2]) << '\n'

        << "g(x, y) = " << g(original_matrix[1][0], x,
            original_matrix[1][1], y,
            original_matrix[1][2]) << '\n';
}
*/