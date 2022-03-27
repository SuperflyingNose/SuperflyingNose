using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab5Math
{
    internal class Program
    {
        static void Main(string[] args)
        {
            int bestI = 0;
            int bestJ = 0;
            double bestPogr = 0;
            double[] minPogr = new double[50];
            for(int i = 1; i < 6; i++)
            {
                for(int j = 1; j < 2; j++)
                {
                    int n = 30 + 10 * j;
                    Console.WriteLine($"Полином {Math.Pow(2, i)} степени\nДля {n} точек");
                    double[] xk = uzli(-3, 3, n);
                    double[] fk = new double[n];
                    for (int k = 0; k < n; k++)
                    {
                        fk[k] = Math.Log(1 + Math.Abs(xk[k]));
                    }
                    double[] P_x = new double[n];
                    double[] a = poiskA((int)Math.Pow(2, i), n, xk, fk, P_x);
                    Console.WriteLine($"-----------------------------------------\n|абсолютная погрешность  {Program.find_delt(fk, P_x).Max()}|\n-----------------------------------------");
                    minPogr[j - 1 + 10 * (i - 1)] = Program.find_delt(fk, P_x).Max();
                    if (minPogr[j - 1 + 10 * (i - 1)] == minPogr.Take(j + 10 * (i - 1)).Min())
                    {
                        bestI = i;
                        bestJ = j;
                        bestPogr = minPogr[j - 1 + 10 * (i - 1)];
                    }
                }
            }

            //Console.WriteLine($"Минимальная абсолютная погрешность = {bestPogr}\n" +
              //  $"Наилучшим полиномом является полином {Math.Pow(2, bestI)} степени для { 10 * bestJ} точек");
            Console.Read();
        }
        static double[] uzli(double bottom, double top, int n)
        {
            double[] mas = new double[n];
            double step = (top - bottom) / (n - 1);
            for (int i = 0; i < n; i++)
            {
                mas[i] = bottom + step * i;
            }
            return mas;
        }
        static double[] poiskA(int m, int n, double[] xk, double[]fk, double[] P_x) // m - степень полинома. n - количесво точек. xk - массив точек. P_x пустой массив, в него вписываются значения аппроксимирующей функции. 
        {
            double[][] mas = new double[m + 1][];
            double[] f = new double[m + 1];
            for (int i = 0; i < m + 1; i++)
            {
                mas[i] = new double[m + 1];
                for (int j = 0; j < m + 1; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        mas[i][j] += Math.Pow(xk[k], i + j);
                    }
                }
                for (int k = 0; k < n; k++)
                {
                    f[i] += fk[k] * Math.Pow(xk[k], i); // fk[x] - значения функции в точке 
                }//Math.Abs(xk[k])    Math.Log(1 + Math.Abs(xk[k]))
            }
            double[] a = gauss(mas, f);
            for(int i = 0; i < xk.Length; i++)
            {
                for(int j = 0; j < m + 1; j++)
                {
                    P_x[i] += a[j] * Math.Pow(xk[i], j);
                }
            }

            //Вывод функции полинома
            /*Console.Write($"P{m}(x) = ");
            for (int i = 0; i < a.Length; i++)
            {
                Console.Write($"{Math.Round(a[i], 4)} * x^{i}");
                if (i != a.Length - 1) { Console.Write(" + "); }
                else { Console.WriteLine(); }
            }*/

            return a;
        }
        static double[] gauss(double[][] mas, double[] F)
        {
            int k = 0;
            while (k < F.Length)
            {
                double o = mas[k][k];
                for (int i = k; i < F.Length; i++)
                {
                    mas[k][i] /= o;
                }
                F[k] /= o;
                for (int i = 1 + k; i < F.Length; i++)
                {
                    o = mas[i][k];
                    F[i] -= F[k] * mas[i][k];
                    for (int j = k; j < mas.Length; j++)
                    {
                        mas[i][j] -= mas[k][j] * o;
                    }
                }
                k++;
            }
            //vivod(mas, F);
            double[] x = new double[F.Length];
            x[F.Length - 1] = F[F.Length - 1];
            double q = F.Length - 2;
            for (int i = F.Length - 2; i >= 0; i--)
            {
                double a = F[i];
                for (int j = F.Length - 1; j > q; j--)
                {
                    a -= mas[i][j] * x[j];
                }
                x[i] = a;
                q--;
            }
            return x;
        }
        static double[] find_delt(double[] reservedF, double[] changedF)
        {
            double[] deltmas = new double[reservedF.Length];
            for (int i = 0; i < deltmas.Length; i++)
            {
                deltmas[i] = Math.Abs(reservedF[i] - changedF[i]);
            }
            return deltmas;
        }
    }

}
