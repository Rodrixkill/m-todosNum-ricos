#include "Diferential.h"
#include <sstream>
#include <iomanip>

namespace met {

	double puntoMedio(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
			y += h * k2;
		}
		return y;
	}

	double eulerModificado(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + h, y + h * k1);
			y += (h / 2)*(k1 + k2);
		}
		return y;
	}

	double heun(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k1);
			y += (h / 4)*(k1 + 3 * k2);
		}
		return y;
	}

	double kutta(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
			double k3 = f.evaluate(x + h, y + h * (-k1 + 2 * k2));
			y += (h / 6)*(k1 + 4 * k2 + k3);
		}
		return y;
	}

	double nystrom(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k1);
			double k3 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k2);
			y += (h / 8)*(2 * k1 + 3 * k2 + 3 * k3);
		}
		return y;
	}

	double rungeKutta(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			double k1 = f.evaluate(x, y);
			double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
			double k3 = f.evaluate(x + h / 2, y + (h / 2 * k2));
			double k4 = f.evaluate(x + h, y + h * k3);
			y += (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
		}
		return y;
	}

	double euler(Parser f, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			y += h * f.evaluate(x, y);
		}
		return y;
	}

	double euler2(Parser f, Parser fx, Parser fy, double x0, double y0, double xf, double h) {
		double y = y0;
		for (double x = x0; x < xf; x += h) {
			y += h * f.evaluate(x, y) + ((h*h) / 2)*(fx.evaluate(x, y) + fy.evaluate(x, y)*f.evaluate(x, y));
		}
		return y;
	}

	string mulitpleValues(Parser f, double x0, double y0, double xf, double h, int func) {
		stringstream ss;

		ss << setprecision(15) << "y( " << x0 << " ) = " << y0 << "\n";

		double y = y0;
		double x = x0;
		double n = (xf - x0) / h;

		for (int i = 0; i < n; i++) {
			if (func == 1) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
				y += h * k2;
			}
			else if (func == 2) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + h, y + h * k1);
				y += (h / 2)*(k1 + k2);;
			}
			else if (func == 3) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k1);
				y += (h / 4)*(k1 + 3 * k2);
			}
			else if (func == 4) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
				double k3 = f.evaluate(x + h, y + h * (-k1 + 2 * k2));
				y += (h / 6)*(k1 + 4 * k2 + k3);;
			}
			else if (func == 5) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k1);
				double k3 = f.evaluate(x + (2 * h) / 3, y + ((2 * h) / 3)*k2);
				y += (h / 8)*(2 * k1 + 3 * k2 + 3 * k3);
			}
			else if (func == 6) {
				double k1 = f.evaluate(x, y);
				double k2 = f.evaluate(x + h / 2, y + (h / 2)*k1);
				double k3 = f.evaluate(x + h / 2, y + (h / 2 * k2));
				double k4 = f.evaluate(x + h, y + h * k3);
				y += (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
			}
			else if (func == 7) {
				y += h * f.evaluate(x, y);
			}
			x += h;

			ss << "y( " << x << " ) = " << y << "\n";
		}

		return ss.str();
	}

	string mulitpleEuler2(Parser f, Parser fx, Parser fy, double x0, double y0, double xf, double h) {
			stringstream ss;

			ss << setprecision(15) << "y( " << x0 << " ) = " << y0 << "\n";

			double y = y0;
			double x = x0;
			double n = (xf - x0) / h;

			for (int i = 0; i < n; i++) {

				y += h * f.evaluate(x, y) + ((h*h) / 2)*(fx.evaluate(x, y) + fy.evaluate(x, y)*f.evaluate(x, y));
				
				x += h;

				ss << "y( " << x << " ) = " << y << "\n";
			}

			return ss.str();
		
	}

	string multSystems(Parser fx, Parser fy, double t0, double x0, double y0, double tf, double h) {
		stringstream ss;
		ss << setprecision(15) << "y( " << t0 << " ) = " << y0 << "          -          " << "x( " << t0 << ") = " << x0 << "\n";

		double y = y0;
		double x = x0;
		double t = t0;
		double n = (tf - t0) / h;

		for (int i = 0; i < n; i++) {

			double auxY = h * fy.evaluate2(t, x, y);

			x += h * fx.evaluate2(t, x, y);
			y += auxY;

			t += h;

			ss << "y( " << t << " ) = " << y << "          -          " << "x( " << t << ") = " << x << "\n";
		}

		return ss.str();
	}

	string systems(Parser fx, Parser fy, double t0, double x0, double y0, double tf, double h) {

		double y = y0;
		double x = x0;
		double t = t0;
		double n = (tf - t0) / h;

		for (int i = 0; i < n; i++) {

			double auxY = h * fy.evaluate2(t, x, y);

			x += h * fx.evaluate2(t, x, y);
			y += auxY;

			t += h;
		}
		stringstream ss;

		ss << setprecision(15) << "y( " << t << " ) = " << y << "\n" << "x( " << t << ") = " << x << "\n";
		return ss.str();
	}



}