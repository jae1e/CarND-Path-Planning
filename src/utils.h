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
namespace Utils
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

///*
//* prediction with constant jerk value,
//* output: { accel, speed, pos }
//*/
//vector<double> matchSpeedPrediction(double jerk, double accel, double speed, double pos, 
//									double max_jerk, double max_accel, double target_speed)
//{
//
//	// max jerk speed
//	vector<double> max_jerk_pred = Utils::constJerkPrediction(ps.max_jerk, ps.prev_dds,
//		car_speed, car_s, pred_dt);
//	double max_jerk_s = max_jerk_pred[2];
//	double max_jerk_ds = max_jerk_pred[1];
//	double max_jerk_dds = max_jerk_pred[0];
//
//	???
//
//	vector<double> traj(3);
//	traj[0] = accel + jerk * dt;
//	traj[1] = speed + accel * dt + 0.5 * jerk * dt * dt;
//	traj[2] = pos + speed * dt + 0.5 * accel * dt * dt + 1.0 / 6.0 * jerk * dt * dt * dt;
//
//	return traj;
//}
//
//double constJerkSpeedDelta(double jerk, double dt)
//{
//	return 0.5 * jerk * dt * dt;
//}
//
///*
//* prediction with constant jerk value, 
//* output: { accel, speed, pos }
//*/
//vector<double> constJerkPrediction(double jerk, double accel, double speed, double pos, double dt)
//{
//	vector<double> traj(3);
//	traj[0] = accel + jerk * dt;
//	traj[1] = speed + accel * dt + 0.5 * jerk * dt * dt;
//	traj[2] = pos + speed * dt + 0.5 * accel * dt * dt + 1.0 / 6.0 * jerk * dt * dt * dt;
//
//	return traj;
//}
//
//vector<double> constJerkTrajectory(double jerk, vector<double> cur_info, int num_trajectory, double move_dt)
//{
//	vector<double> traj(num_trajectory);
//
//	double pos = cur_info[0];
//	double speed = cur_info[1];
//	double accel = cur_info[2];
//
//	// fill trajectory
//	for (int it = 0; it < num_trajectory; ++it)
//	{
//		vector<double> pred = constJerkPrediction(jerk, accel, speed, pos, (it + 1) * move_dt);
//		traj[it] = pred[2];
//	}
//
//	return traj;
//}

double calculateAlpha(double accel, double max_accel, double coeff1, double coeff2)
{
	return 0.5 * (1.0 + coeff1 * exp(-coeff2 * accel / max_accel));
}

double efficientDuration(vector< double> start, vector <double> end, double alpha)
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
	*/
	double eps = 1e-2;
	return (end[0] - start[0]) / (start[1] + alpha * (end[1] - start[1]) + eps);
}

double efficientDistance(double start_pos, double start_speed, double target_speed, double T, double alpha)
{
	/*
	Calculate the Jerk Minimizing Trajectory that connects the initial state
	to the final state in time T.

	INPUTS

	start - the vehicles start location given as a length three array
	corresponding to initial values of [s, s_dot]

	end   - the desired end state for vehicle. Like "start" this is a
	length three array.

	T     - The duration, in seconds, over which this maneuver should occur.
	*/
	
	return start_pos + (start_speed + 1 * (target_speed - start_speed)) * T;
	return start_pos + (start_speed + alpha * (target_speed - start_speed)) * T;
}

vector<double> jerkMinimizingCoeffs(vector< double> start, vector <double> end, double T)
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

vector<double> jerkMinimizingTrajectory(vector<double> cur_info, vector<double> tgt_info, 
										int num_trajectory, double pred_dt, double iter_dt)
{
	vector<double> traj(num_trajectory);

	// polynomial fitting
	auto coeffs = jerkMinimizingCoeffs(cur_info, tgt_info, pred_dt);

	// fill trajectory
	for (int it = 0; it < num_trajectory; ++it)
	{
		double val = 0;
		for (int ic = 0; ic < coeffs.size(); ++ic)
		{
			val += coeffs[ic] * pow((it + 1) * iter_dt, ic);
		}

		traj[it] = val;
	}

	return traj;
}

/*
 * prec_s_dist: s position of preceding car - s position of ego car
 * prec_ds: speed of preceding car
 * delta_t: timespan of interest
 */
double prececedingDistPred(double prec_s_dist, double prec_ds, double delta_t)
{
	return prec_s_dist + prec_ds * delta_t;
}

} // namespace Utils
	
}

#endif /* utils_h */
