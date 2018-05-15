#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define REPEATS 100.0

// Варианты критической области
#define DOUBLE_SIDED 1
#define ONE_SIDED_L 2
#define ONE_SIDED_R 3

void W(float *x, int n, float *y, int m, float alpha, int critical_area);	// Критерий Вилкоксона
void Q(float *x, int m, float *y, int n, float alpha, float critical_value);// Критерий Крамера-Уэлча
void t(float *x, int m, float *y, int n, float alpha, float critical_value);// Критерий Стьюдента

int main()
{
	// Задание выборок
	int m = 10, n = 10;
	float x[] = {76, 49, 83, 72, 64, 73, 45, 52, 26, 66};
	float y[] = {46, 71, 78, 58, 96, 77, 85, 48, 69, 99};

	// Тестовый пример
//	int m = 4, n = 6;
//	float x[] = {1, 2, 3, 4};
//	float y[] = {1, 2, 3, 4, 5, 6};
	W(x, m, y, n, 0.05f, DOUBLE_SIDED);
	Q(x, m, y, n, 0.05f, 1.96f);
	t(x, m, y, n, 0.05f, 10.0f);

	return 0;
}

void W(float *x, int m, float *y, int n, float alpha, int critical_area)
{
	int *r;			// Массив с текущими рангами xi
	float *z;		// Объединенный массив
	int sum, min_sum, max_sum;
	int i, j;
	float p;		// Вероятность отклонения нулевой гипотезы
	int *ns;		// Массив с частотами появления сумм рангов
	
	for(min_sum = 0, i = 1; i <= m; i++)
		min_sum += i;

	for(max_sum = 0, i = n + 1; i <= m + n; i++)
		max_sum += i;

	r = new int[m];
	ns = new int[max_sum - min_sum + 1];
	z = new float[m + n];
	
	for(i = 0; i < m; i++)
		r[i] = i + 1;

	for(i = 0; i <= max_sum - min_sum; i++)
		ns[i] = 0;
	
	// Построение закона распределения для суммы рангов
	i = m - 1;
	while(1)
	{
		for(sum = 0, j = 0; j < m; j++)
			sum += r[j];

		ns[sum - min_sum] += 1;

		if(r[m - 1] == m + n)
		{
			i = m - 2;
			while(r[i + 1] - r[i] < 2 && i != -1)
				i--;
			if(i == -1) break;

			r[i]++;

			for(j = i + 1; j < m; j++)
				r[j] = r[j - 1] + 1;
		}
		else
			r[m - 1]++;
	}

	// Общее количество комбинаций
	for(i = 0; i <= max_sum - min_sum; i++)
		j += ns[i];
		
	// Вычисление доверительного интервала
	p = 0;
	
	if(critical_area == DOUBLE_SIDED)
	{
		for(i = 0; p < alpha; i++)
			p += (float) (ns[i] + ns[max_sum - min_sum - i]) / j;
		i--;
		printf("Interval: [%d; %d]\n", min_sum + i, max_sum - i);
	}

	if(critical_area == ONE_SIDED_L)
	{
		for(i = 0; p < alpha; i++)
			p += (float) ns[max_sum - min_sum - i] / j;
		i--;
		printf("Interval: [-infinity; %d]\n", max_sum - i);
	}
	
	if(critical_area == ONE_SIDED_R)
	{
		for(i = 0; p < alpha; i++)
			p += (float) ns[i] / j;
		i--;
		printf("Interval: [%d; +infinity]\n", min_sum + i);
	}

	printf("M(S) = %.1f\n", (float) m * (m + n + 1) / 2);


	// Нахождение наблюдаемого значения критерия
	srand((unsigned)time(NULL));

	
	sum = 0;
	for(int k = 0; k < REPEATS; k++)
	{
        for(i = 0; i < m; i++)
			z[i] = x[i];
	
		for(i = 0; i < n; i++)
			z[i + m] = y[i];
	
		for(i = 0; i < m; i++)
			r[i] = i + 1;

		for(i = 0; i < m + n; i++)
			for(j = 0; j < m + n - 1; j++)
				if(z[j] > z[j + 1] || (z[j] == z[j + 1] && rand()%2 == 0))
				{
					int k1, k2;
				
					for(k1 = 0; k1 < m; k1++)
						if(r[k1] == j + 1) break;
				
					for(k2 = 0; k2 < m; k2++)
						if(r[k2] == j + 2) break;
				
					if(k1 < m) r[k1]++;
					if(k2 < m) r[k2]--;
					
					p = z[j];
					z[j] = z[j + 1];
					z[j + 1] = z[j];
				}
					
		for(i = 0; i < m; i++)
			sum += r[i];
	}
	
	printf("S = %f", sum / REPEATS);
	
	getch();

	return;
}

void Q(float *x, int m, float *y, int n, float alpha, float critical_value)
{
	float z_av = 0, s = 0, Q;
	int i;

	for(i = 0; i < m; i++)
		z_av += x[i] - y[i];

	z_av /= m;

	for(i = 0; i < m; i++)
		s += (x[i] - y[i] - z_av) * (x[i] - y[i] - z_av);

	s /= m - 1;

	s = powf(s, 0.5f);

	Q = powf((float) n, 0.5f) * z_av / s;

	printf("\nQ = %f, Qcr = %f", Q, critical_value);
	getch();
	return;
}

void t(float *x, int m, float *y, int n, float alpha, float critical_value)
{
	int i;
	float x_av = 0, y_av = 0, s2x = 0, s2y = 0, t;
	
	for(i = 0; i < m; i++)
		x_av += x[i];
	x_av /= m;

	for(i = 0; i < m; i++)
		s2x += (x[i] - x_av) * (x[i] - x_av);
	s2x /= m - 1;

	for(i = 0; i < n; i++)
		y_av += y[i];
	y_av /= n;

	for(i = 0; i < n; i++)
		s2y += (y[i] - y_av) * (y[i] - y_av);
	s2y /= n - 1;

	t = (x_av - y_av) / powf((m - 1) * s2x + (n - 1) * s2y, 0.5f) *
		powf((float) m * n * (m + n - 2) / (m + n), 0.5f);

	printf("\nt = %f, tcr = %f", t, critical_value);
	getch();
	return;
}