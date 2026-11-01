// Linear Equation Solver: Gauss elimination, Gauss-Jordan and Cramer's rule
// Libraries: iostream (input/output), cmath (fabs), vector (matrix of any size)

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

// a matrix is a vector of rows, and each row is a vector of numbers
typedef vector< vector<double> > Matrix;

// numbers smaller than this are treated as zero
const double EPS = 0.000000001;

// ===== 1. INPUT ====================================================

// reads a whole number that is at least "low", asks again if invalid
int readInt(const char prompt[], int low) {
    int v;
    while (true) {
        cout << prompt;
        if (cin >> v && v >= low) return v;
        cout << "Invalid, enter again\n";
        cin.clear();
        cin.ignore(10000, '\n');
    }
}

// reads y or n, returns true for y
bool readYesNo(const char prompt[]) {
    char c;
    while (true) {
        cout << prompt;
        cin >> c;
        if (c == 'y' || c == 'Y') return true;
        if (c == 'n' || c == 'N') return false;
        cout << "Invalid, enter y or n\n";
    }
}

// Function 1: take input and build the augmented matrix (m rows, n+1 columns)
void inputSystem(Matrix &A, int &m, int &n) {
    m = readInt("Enter the number of equations (m): ", 1);
    n = readInt("Enter the number of unknowns (n): ", 1);

    A.assign(m, vector<double>(n + 1, 0.0));

    cout << "\nEnter each equation as " << n
         << " coefficients followed by the constant.\n";
    cout << "Example: 2 1 -1 8  means  2x + y - z = 8\n\n";

    for (int i = 0; i < m; i++) {
        bool ok = false;
        while (!ok) {
            cout << "Equation " << i + 1 << ": ";
            ok = true;
            for (int j = 0; j <= n; j++) {
                if (!(cin >> A[i][j])) {
                    cin.clear();
                    cin.ignore(10000, '\n');
                    cout << "Invalid, enter the whole row again\n";
                    ok = false;
                    break;
                }
            }
        }
    }
}

// ===== 2. CORE MATH (shared by all methods) ========================

// Function 2: row echelon form (any size). Returns the rank; fills pivotCol and swaps.
int rowEchelon(Matrix &A, int m, int n, vector<int> &pivotCol, int &swaps) {
    swaps = 0;
    pivotCol.clear();
    int row = 0;

    for (int col = 0; col < n && row < m; col++) {
        // find the row (from "row" downward) with the largest number in this column
        int best = row;
        for (int i = row + 1; i < m; i++) {
            if (fabs(A[i][col]) > fabs(A[best][col])) best = i;
        }

        // no pivot in this column, go to the next column
        if (fabs(A[best][col]) < EPS) continue;

        if (best != row) {
            for (int j = 0; j <= n; j++) {
                double temp = A[best][j];
                A[best][j] = A[row][j];
                A[row][j] = temp;
            }
            swaps++;
        }

        for (int i = row + 1; i < m; i++) {
            double factor = A[i][col] / A[row][col];
            for (int j = col; j <= n; j++) {
                A[i][j] = A[i][j] - factor * A[row][j];
            }
            A[i][col] = 0;
        }

        pivotCol.push_back(col);
        row++;
    }

    return row;
}

// returns 0 = unique solution, 1 = no solution, 2 = infinite solutions
int checkResult(const Matrix &A, int m, int n, int rank) {
    // a row like 0 0 0 | 5 means no solution
    for (int i = rank; i < m; i++) {
        if (fabs(A[i][n]) > EPS) return 1;
    }
    if (rank < n) return 2;
    return 0;
}

// continues after row echelon: makes every pivot 1 and makes zeros above it
void reduceFurther(Matrix &A, int n, int rank, const vector<int> &pivotCol) {
    for (int k = rank - 1; k >= 0; k--) {
        int c = pivotCol[k];

        double pivot = A[k][c];
        for (int j = 0; j <= n; j++) {
            A[k][j] = A[k][j] / pivot;
        }

        for (int i = 0; i < k; i++) {
            double factor = A[i][c];
            for (int j = 0; j <= n; j++) {
                A[i][j] = A[i][j] - factor * A[k][j];
            }
        }
    }
}

// Function 3: determinant of an n x n matrix (works on a copy, uses row echelon).
double determinant(Matrix B, int n) {
    vector<int> pivotCol;
    int swaps;
    int rank = rowEchelon(B, n, n, pivotCol, swaps);
    if (rank < n) return 0;

    double det = 1;
    for (int i = 0; i < n; i++) det = det * B[i][i];
    if (swaps % 2 == 1) det = -det;  // each row swap flips the sign
    return det;
}

// ===== 3. PRINTING THE ANSWER ======================================

void printSolution(const vector<double> &x) {
    cout << "\nSolution:\n";
    for (int i = 0; i < (int)x.size(); i++) {
        double v = x[i];
        if (fabs(v) < 0.00005) v = 0;  // avoid printing -0
        cout << "x" << i + 1 << " = " << v << "\n";
    }
}

// prints infinite solutions with free variables t1, t2, ... (matrix must be reduced)
void printInfinite(const Matrix &A, int n, int rank, const vector<int> &pivotCol) {
    // rowOfCol[c] = the pivot row of column c, or -1 if column c has no pivot
    vector<int> rowOfCol(n, -1);
    for (int k = 0; k < rank; k++) rowOfCol[pivotCol[k]] = k;

    vector<int> tNumber(n, 0);
    int freeCount = 0;
    for (int c = 0; c < n; c++) {
        if (rowOfCol[c] == -1) {
            freeCount++;
            tNumber[c] = freeCount;
        }
    }

    cout << "\nInfinite solutions (" << freeCount
         << " free variable(s), t can be any number):\n";

    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = ";

        if (rowOfCol[i] == -1) {
            cout << "t" << tNumber[i] << "\n";   // a free unknown
        } else {
            int k = rowOfCol[i];
            double constant = A[k][n];
            if (fabs(constant) < 0.00005) constant = 0;
            cout << constant;

            for (int j = 0; j < n; j++) {
                if (rowOfCol[j] == -1 && fabs(A[k][j]) > EPS) {
                    double coef = A[k][j];
                    if (coef > 0) cout << " - " << coef << "*t" << tNumber[j];
                    else          cout << " + " << -coef << "*t" << tNumber[j];
                }
            }
            cout << "\n";
        }
    }
}

// ===== 4. THE THREE METHODS ========================================

// Method 1: Gauss elimination
void gaussElimination(Matrix &A, int m, int n) {
    vector<int> pivotCol;
    int swaps;
    int rank = rowEchelon(A, m, n, pivotCol, swaps);
    int outcome = checkResult(A, m, n, rank);

    if (outcome == 1) {
        cout << "\nNo solution (the equations contradict each other).\n";
        return;
    }
    if (outcome == 2) {
        reduceFurther(A, n, rank, pivotCol);
        printInfinite(A, n, rank, pivotCol);
        return;
    }

    // back substitution: a unique solution means row i has its pivot in column i
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = A[i][n];
        for (int j = i + 1; j < n; j++) {
            sum = sum - A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
    printSolution(x);
}

// Method 2: Gauss-Jordan
void gaussJordan(Matrix &A, int m, int n) {
    vector<int> pivotCol;
    int swaps;
    int rank = rowEchelon(A, m, n, pivotCol, swaps);
    int outcome = checkResult(A, m, n, rank);

    if (outcome == 1) {
        cout << "\nNo solution (the equations contradict each other).\n";
        return;
    }

    reduceFurther(A, n, rank, pivotCol);

    if (outcome == 2) {
        printInfinite(A, n, rank, pivotCol);
        return;
    }

    // the left side is now the identity, so the answers are in the last column
    vector<double> x(n, 0.0);
    for (int i = 0; i < n; i++) x[i] = A[i][n];
    printSolution(x);
}

// Method 3: Cramer's rule (square systems only)
void cramersRule(const Matrix &A, int m, int n) {
    if (m != n) {
        cout << "\nCramer's rule needs a square system (m = n). "
                "Cannot use it here.\n";
        return;
    }

    double det = determinant(A, n);
    if (fabs(det) < EPS) {
        cout << "\nDeterminant is 0, so Cramer's rule cannot give a unique solution.\n"
                "Try Gauss elimination or Gauss-Jordan, which can tell you whether\n"
                "there is no solution or infinite solutions.\n";
        return;
    }

    vector<double> x(n, 0.0);
    for (int j = 0; j < n; j++) {
        Matrix D = A;
        for (int i = 0; i < n; i++) D[i][j] = A[i][n];
        x[j] = determinant(D, n) / det;
    }
    printSolution(x);
}

// ===== 5. MENU =====================================================

int chooseMethod() {
    cout << "\n===== Linear Equation Solver =====\n";
    cout << "1. Gauss elimination\n";
    cout << "2. Gauss-Jordan\n";
    cout << "3. Cramer's rule\n";

    int choice = 0;
    while (true) {
        cout << "Enter your choice (1-3): ";
        if (cin >> choice && choice >= 1 && choice <= 3) return choice;
        cout << "Invalid choice, choose again\n";
        cin.clear();
        cin.ignore(10000, '\n');
    }
}

// ===== 6. MAIN PROGRAM =============================================

int main() {
    bool again = true;
    while (again) {
        int choice = chooseMethod();

        if (choice == 3) {
            cout << "\nLimitation: Cramer's rule works only for square systems "
                    "(equations = unknowns)\nwith a non-zero determinant.\n";
            if (!readYesNo("Do you want to continue? (y/n): ")) {
                continue;
            }
        }

        Matrix A;
        int m, n;
        inputSystem(A, m, n);

        switch (choice) {
            case 1: gaussElimination(A, m, n); break;
            case 2: gaussJordan(A, m, n); break;
            case 3: cramersRule(A, m, n); break;
        }

        again = readYesNo("\nSolve another system? (y/n): ");
    }

    cout << "Goodbye!\n";
    return 0;
}