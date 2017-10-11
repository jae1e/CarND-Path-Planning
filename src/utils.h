//
//  trajectory.h
//  Path_Planning
//
//  Created by Jaeil Park on 2017. 10. 10..
//
//

#ifndef utils_h
#define utils_h

#include <vector>
#include "spline.h"

using namespace std;

namespace
{
namespace Trajectory
{
	
void interpolate_points(vector<double> along_points, vector<double> input_points, int ratio, vector<double>& output_points)
{
	output_points.clear();
	
	tk::spline spline;
	spline.set_points(along_points, input_points);
	
	for (int i = 0; i < along_points.size() - 1; ++i)
	{
		double p1 = along_points[i];
		double p2 = along_points[i+1];
		double dp = (p2 - p1) / ratio;
		
		for (int j = 0; j < ratio; ++j)
		{
			output_points.push_back(spline(p1 + j * dp));
		}
	}
	
	output_points.push_back(along_points.back());
}

} // namespace Trajectory
	
namespace CostFunction
{

/*
 * min_ds: s distance between ego car and preceding car
 * dspeed: speed of preceding car - speed of ego car
 * delta_t: timespan of interest
 */
double lane_cost(double min_ds, double dspeed, double delta_t)
{
	return -(min_ds + dspeed * delta_t);
}

	
} // namespace CostFunction
	
}

#endif /* utils_h */
