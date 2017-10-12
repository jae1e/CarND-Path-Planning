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
#include "Eigen-3.3/Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace
{
namespace Trajectory
{
	
void interpolatePoints(vector<double> along_points, vector<double> input_points, int ratio, vector<double>& output_points)
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

vector<double> JMT(vector< double> start, vector <double> end, double T)
{
	/*
	Calculate the Jerk Minimizing Trajectory that connects the initial state
	to the final state in time T.

	INPUTS

	start - the vehicles start location given as a length three array
	corresponding to initial values of [s, s_dot, s_double_dot]

	end   - the desired end state for vehicle. Like "start" this is a
	length three array.

	T     - The duration, in seconds, over which this maneuver should occur.

	OUTPUT
	an array of length 6, each value corresponding to a coefficent in the polynomial
	s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

	EXAMPLE

	> JMT( [0, 10, 0], [10, 10, 0], 1)
	[0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
	*/

	MatrixXd A = MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
		3 * T*T, 4 * T*T*T, 5 * T*T*T*T,
		6 * T, 12 * T*T, 20 * T*T*T;

	MatrixXd B = MatrixXd(3, 1);
	B << end[0] - (start[0] + start[1] * T + .5*start[2] * T*T),
		end[1] - (start[1] + start[2] * T),
		end[2] - start[2];

	MatrixXd Ai = A.inverse();

	MatrixXd C = Ai*B;

	vector <double> result = { start[0], start[1], .5*start[2] };
	for (int i = 0; i < C.size(); i++)
	{
		result.push_back(C.data()[i]);
	}

	return result;

}

vector<double> generateTrajectory(vector<double> cur_info, vector<double> tgt_info, int num_trajectory, double pred_dt, double move_dt)
{
	// polynomial fitting
	auto coeffs = JMT(cur_info, tgt_info, pred_dt);

	// generate trajectory
	vector<double> traj(num_trajectory);

	for (int it = 0; it < num_trajectory; ++it)
	{
		double val = 0;
		for (int ic = 0; ic < coeffs.size(); ++ic)
		{
			val += coeffs[ic] * pow(move_dt, ic);
		}

		traj[it] = val;
	}

	return traj;
}

} // namespace Trajectory
	
namespace CostFunction
{

/*
 * min_ds: s distance between ego car and preceding car
 * dspeed: speed of preceding car - speed of ego car
 * delta_t: timespan of interest
 */
double laneCost(double min_ds, double dspeed, double delta_t)
{
	return -(min_ds + dspeed * delta_t);
}

	
} // namespace CostFunction
	
}

#endif /* utils_h */
