

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <chrono>
#include <vector>
#include "../../sparse_matrices.hpp"


using namespace std;
using namespace std::chrono;

int main()
{
    srand((size_t)time(NULL));             // initialization of randomizer

    char c = ' ';
    size_t N;
    size_t P;
    while (c != 'q') {
        cout << "q - quit\nt - test\n";
        cin >> c;

        if (c == 't') {
            cout << "N:";
            cin >> N;
            cout << "P/1000:";
            cin >> P;

			high_resolution_clock::time_point t1, t2;
			duration<double> time_span;

			t1 = high_resolution_clock::now();
            input_storage_scheme<double> ISS(N,N);
            for (size_t i = 0; i < N; ++i)
                ISS.add_element(static_cast<double>(rand()), i, i);
			for (size_t row = 0; row < N; ++row)
				for (size_t col = 0; col < N; ++col)
					if (row != col && static_cast<size_t>(rand() % 1000) < P) {
						ISS.add_element(static_cast<double>(rand()), row, col);
					}
			t2 = high_resolution_clock::now();
			time_span = duration_cast<duration<double>>(t2 - t1);
			cout << "ISS init: " << time_span.count() << "s\n";

			t1 = high_resolution_clock::now();
			dynamic_storage_scheme<double> DSS(ISS, 50);
			t2 = high_resolution_clock::now();
			time_span = duration_cast<duration<double>>(t2 - t1);
			cout << "DSS init: " << time_span.count() << "s\n";

			t1 = high_resolution_clock::now();
			try {
				DSS.LU_decomposition(MARKOWITZ_COST, 5, 5, 0.00001);
			}
			catch(exception& e) {
				cout << e.what() << endl;
				continue;
			}
			t2 = high_resolution_clock::now();
			time_span = duration_cast<duration<double>>(t2 - t1);
			cout << "DSS LU_decomposition: " << time_span.count() << "s\n";

			vector<double> x(N), b(N), r(N);
			for (size_t i = 0; i < N; ++i)
                b[i] = static_cast<double>(rand());

			t1 = high_resolution_clock::now();
			DSS.solve_LU(x.data(), b.data());
			DSS.iterative_refinement(ISS, x.data(), b.data(), 0.0000000001, 100);
			t2 = high_resolution_clock::now();
			time_span = duration_cast<duration<double>>(t2 - t1);
			cout << "DSS iterative refinement: " << time_span.count() << "s\n";

			ISS.count_rasidual_vector(x.data(), b.data(), r.data());
			cout << "residual vector norm: " << vector_norm(r.data(), N) << endl;


        }

        system("pause");
    }

    

return 0;
}
