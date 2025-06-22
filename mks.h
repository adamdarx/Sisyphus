/*
    filename: mks.h
    本文件是一个度规定义的示例。任何度规函数(不依赖于格点)的类型必须是MetricComponent。
    必须声明度规名称。
    如果自定义度规需要额外参数，需要在自定义文件中进行明确指定。
    度规名称：
        Modified Kerr Schwarzchild
    额外参数：
        a: 自旋
        h
*/
#pragma once
#include <cmath>
#include "mks.h"

namespace MKS
{
	const char* metric_name = "Modified Kerr Schwarzchild";  // 度规名称
	const double a = (0.9);
	const double h = (0.);
	const double PI = (3.1415926535897);

	void init_metricFunc(MetricComponent metricFunc[4][4])
	{
		// 度规初始化 
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				metricFunc[i][j] = [](double x0, double x1, double x2, double x3) {return 0.; };
		// 自定义度规
		metricFunc[0][0] = [](double x0, double x1, double x2, double x3) {return -1. + 2. * exp(x1) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.)); };
		metricFunc[0][1] = [](double x0, double x1, double x2, double x3) {return (2. * exp(x1) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.))) * exp(x1); };
		metricFunc[0][3] = [](double x0, double x1, double x2, double x3) {return (-2. * a * exp(x1) * pow(sin(x2 + 0.5 * h * sin(2 * x2)), 2.) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.))); };
		metricFunc[1][0] = metricFunc[0][1];
		metricFunc[1][1] = [](double x0, double x1, double x2, double x3) {return (1. + 2. * exp(x1) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.))) * pow(exp(x1), 2); };
		metricFunc[1][3] = [](double x0, double x1, double x2, double x3) {return (-a * pow(sin(x2 + 0.5 * h * sin(2 * x2)), 2.) * (1. + 2. * exp(x1) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.)))) * exp(x1); };
		metricFunc[2][2] = [](double x0, double x1, double x2, double x3) {return (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.)) * pow((1. + h * cos(2. * x2)), 2); };
		metricFunc[3][0] = metricFunc[0][3];
		metricFunc[3][1] = metricFunc[1][3];
		metricFunc[3][3] = [](double x0, double x1, double x2, double x3) {return (pow(sin(x2 + 0.5 * h * sin(2 * x2)), 2.) * ((pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.)) + pow(a * sin(x2 + 0.5 * h * sin(2 * x2)), 2.) * (1. + 2. * exp(x1) / (pow(exp(x1), 2.) + pow(a * cos(x2 + 0.5 * h * sin(2 * x2)), 2.))))); };
	}
}