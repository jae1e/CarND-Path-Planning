#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "utils.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

enum class BehaviorStatus
{
	KeepLane = 0,
	ChangeLane = 1
};

struct ParameterSet
{
public: // general parameters
	double pred_dt = 0.2;
	
	double move_dt = 0.02;

	int num_pred = 10;
	
	double lane_width = 4.0;
	
	double max_speed = 22; // 50 MPH = 22.352 m/s;
	
	double max_accel = 10.0;
	
	double max_jerk = 10.0;
	
	double safety_distance = 50.0;
	
	double lane_change_distance = 70.0;
	
	double lane_change_cost_coeff = 1.2;

	double lane_change_dt = 2.0;

	double lane_change_speed_decay = 0.2;
	
	int num_lanes = 3;

public: // interpoaltion parameters
	int num_src_waypoints = 10;
	
	int ratio_interpolation = 10;
	
public: // status parameters
	int current_lane_id = 1;
	
	int target_lane_id = 1;
	
	BehaviorStatus current_status = BehaviorStatus::KeepLane;
};

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
	
	// parameter set
	ParameterSet ps;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &ps](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			
			json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	//////////////////////////////////////////////////////////////////////////////////////////////////////
			// interpolate waypoints
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			vector<double> src_waypoints_s, src_waypoints_x, src_waypoints_y, src_waypoints_dx, src_waypoints_dy;
			vector<double> waypoints_s, waypoints_x, waypoints_y, waypoints_dx, waypoints_dy;
			
			{
				int next_waypoint_id = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
				
				int begin_waypoint_id = next_waypoint_id - ps.num_src_waypoints;
				int end_waypoint_id = next_waypoint_id + ps.num_src_waypoints;
				if (begin_waypoint_id < 0)
				{
					begin_waypoint_id = 0;
					end_waypoint_id = begin_waypoint_id + ps.num_src_waypoints;
				}
				else if (end_waypoint_id > (int)map_waypoints_x.size())
				{
					end_waypoint_id = (int)map_waypoints_x.size();
					begin_waypoint_id = end_waypoint_id - ps.num_src_waypoints;
				}
				
				for (int wid = begin_waypoint_id; wid < end_waypoint_id; ++wid)
				{
					src_waypoints_s.push_back(map_waypoints_s[wid]);
					src_waypoints_x.push_back(map_waypoints_x[wid]);
					src_waypoints_y.push_back(map_waypoints_y[wid]);
					src_waypoints_dx.push_back(map_waypoints_dx[wid]);
					src_waypoints_dy.push_back(map_waypoints_dy[wid]);
				}
				
				for (int i = 0; i < src_waypoints_s.size() - 1; ++i)
				{
					double p1 = src_waypoints_s[i];
					double p2 = src_waypoints_s[i+1];
					double dp = (p2 - p1) / ps.ratio_interpolation;
					
					for (int j = 0; j < ps.ratio_interpolation; ++j)
					{
						waypoints_s.push_back(p1 + j * dp);
					}
				}
				waypoints_s.push_back(src_waypoints_s.back());
				Trajectory::interpolatePoints(src_waypoints_s, src_waypoints_x, ps.ratio_interpolation, waypoints_x);
				Trajectory::interpolatePoints(src_waypoints_s, src_waypoints_y, ps.ratio_interpolation, waypoints_y);
				Trajectory::interpolatePoints(src_waypoints_s, src_waypoints_dx, ps.ratio_interpolation, waypoints_dx);
				Trajectory::interpolatePoints(src_waypoints_s, src_waypoints_dy, ps.ratio_interpolation, waypoints_dy);
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// make decision
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			// minimum distance in s of each lane
			vector<double> lane_min_ds(ps.num_lanes, 1000);
			vector<double> lane_dspeed(ps.num_lanes, 0);

			{
				// update current lane and current status
				ps.current_lane_id = int(car_d / ps.lane_width);
				
				if (ps.current_status == BehaviorStatus::ChangeLane)
				{
					if (ps.target_lane_id == ps.current_lane_id)
					{
						ps.current_status = BehaviorStatus::KeepLane;
					}
				}
				
				// consider lane change only when the current status is keep lane and d is within current lane's safe range
				double current_lane_d = (0.5 + ps.current_lane_id) * ps.lane_width;
				double lane_safe_d_width = 0.5 * ps.lane_width;
				
				if (ps.current_status == BehaviorStatus::KeepLane
					&& car_d > current_lane_d - 0.5 * lane_safe_d_width && car_d < current_lane_d + 0.5 * lane_safe_d_width)
				{
					// fill ds and dspeed of each lane
					for (auto sfit = sensor_fusion.begin(); sfit != sensor_fusion.end(); ++sfit)
					{
						// int id = sfit->at(0);
						// double x = sfit->at(1), y = sfit->at(2);
						double vx = sfit->at(3), vy = sfit->at(4);
						double s = sfit->at(5), d = sfit->at(6);
						
						int lane_id = (int)(d / ps.lane_width);
						if (lane_id == ps.current_lane_id)
						{
							if (s > car_s && s - car_s < lane_min_ds[lane_id])
							{
								lane_min_ds[lane_id] = s - car_s;
								
								double speed = sqrt(vx * vx + vy * vy);
								lane_dspeed[lane_id] = speed - car_speed;
							}
						}
					}
					
					// if current lane's min distance and delta speed is not enough,
					if (lane_min_ds[ps.current_lane_id] < ps.lane_change_distance
						&& lane_dspeed[ps.current_lane_id] < 0)
					{
						// calculate lane cost
						double current_cost = CostFunction::laneCost(lane_min_ds[ps.current_lane_id],
																	  lane_dspeed[ps.current_lane_id],
																	  ps.pred_dt);
						
						double left_cost = 0, right_cost = 0;
						if (ps.current_lane_id > 0)
						{
							left_cost = CostFunction::laneCost(lane_min_ds[ps.current_lane_id - 1],
																lane_dspeed[ps.current_lane_id - 1],
																ps.pred_dt);
						}
						if (ps.current_lane_id < ps.num_lanes - 1)
						{
							right_cost = CostFunction::laneCost(lane_min_ds[ps.current_lane_id + 1],
																 lane_dspeed[ps.current_lane_id + 1],
																ps.pred_dt);
						}
						
						// make decision
						int target_lane_id = right_cost < left_cost ? ps.current_lane_id + 1 : ps.current_lane_id - 1;
						double min_cost = min(right_cost, left_cost);
						if (current_cost * ps.lane_change_cost_coeff > min_cost)
						{
							ps.target_lane_id = target_lane_id;
							ps.current_status = BehaviorStatus::ChangeLane;
						}
					}
				}
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// target points generation
			//////////////////////////////////////////////////////////////////////////////////////////////////////

			double target_s = 0, target_ds = 0, target_dds = 0;
			double target_d = 0, target_dd = 0, target_ddd = 0;

			{
				double lane_change_speed = car_speed * (1.0 - ps.lane_change_speed_decay);

				// check viability : no cars within safety range while changing lane
				if (ps.current_status == BehaviorStatus::ChangeLane)
				{
					bool lane_change_viable = true;

					double total_dt = ps.lane_change_dt;
					double target_s = car_s + lane_change_speed * total_dt;

					double min_safe_s = target_s - ps.safety_distance;
					double max_safe_s = target_s + ps.safety_distance;

					for (auto sfit = sensor_fusion.begin(); sfit != sensor_fusion.end(); ++sfit)
					{
						// int id = sfit->at(0);
						// double x = sfit->at(1), y = sfit->at(2);
						double vx = sfit->at(3), vy = sfit->at(4);
						double s = sfit->at(5), d = sfit->at(6);
						int lane_id = int(d / ps.lane_width);
						if (lane_id != ps.current_lane_id || lane_id != ps.target_lane_id)
						{
							continue;
						}

						double speed = sqrt(vx * vx + vy * vy);
						double end_s = s + speed * total_dt;

						if (end_s > min_safe_s && end_s < max_safe_s)
						{
							lane_change_viable = false;
							break;
						}
					}

					if (!lane_change_viable)
					{
						ps.current_status = BehaviorStatus::KeepLane;
						ps.target_lane_id = ps.current_lane_id;
					}
				}

				// set target
				if (ps.current_status == BehaviorStatus::ChangeLane)
				{
					target_s = car_s + lane_change_speed * ps.pred_dt;
					target_ds = lane_change_speed;
				}
				else // if (ps.current_status == BehaviorStatus::KeepLane)
				{
					double predicted_preceding_car_s = lane_min_ds[ps.current_lane_id] + lane_dspeed[ps.current_lane_id] * ps.pred_dt;
					target_s = predicted_preceding_car_s - ps.safety_distance;
					double keep_dist_speed = (target_s - car_s) / ps.pred_dt;
					target_ds = min(ps.max_speed, keep_dist_speed);
				}
				target_dds = 0;

				target_d = (0.5 + ps.target_lane_id) * ps.lane_width;
				target_dd = 0;
				target_ddd = 0;
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// generate trajectory
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			{
				double dt = ps.move_dt;
				double s_cur = 0, ds_cur = 0, dds_cur = 0;
				double d_cur = 0, dd_cur = 0, ddd_cur = 0;

				int num_prev_path = previous_path_x.size();

				// calculate current s and d info
				if (previous_path_x.size() < 5)
				{
					s_cur = car_s;
					ds_cur = car_speed;
					dds_cur = 0;

					d_cur = car_d;
					dd_cur = 0;
					ddd_cur = 0;
				}
				else
				{
					vector<double> x_prev(3);
					vector<double> y_prev(3);
					vector<double> s_prev(3);
					vector<double> d_prev(3);

					// get frenet coord of previous path
					for (int i = 0; i < 3; ++i)
					{
						double x1 = previous_path_x[num_prev_path - 3 + i - 1];
						double x2 = previous_path_x[num_prev_path - 3 + i];
						double y1 = previous_path_y[num_prev_path - 3 + i - 1];
						double y2 = previous_path_y[num_prev_path - 3 + i];

						double angle = atan2(y2 - y1, x2 - x1);
						vector<double> frenet = getFrenet(x2, y2, angle, waypoints_x, waypoints_y);

						x_prev[i] = x2;
						y_prev[i] = y2;
						s_prev[i] = frenet[0];
						d_prev[i] = frenet[1];
					}

					s_cur = s_prev[2];
					ds_cur = (s_prev[2] - s_prev[1]) / dt;
					dds_cur = (ds_cur - ((s_prev[1] - s_prev[0]) / dt)) / dt;

					d_cur = d_prev[2];
					dd_cur = (d_prev[2] - d_prev[1]) / dt;
					ddd_cur = (dd_cur - ((d_prev[1] - d_prev[0]) / dt)) / dt;
				}

				vector<double> cur_s_info = { s_cur, ds_cur, dds_cur };
				vector<double> tgt_s_info = { target_s, target_ds, target_dds };
				vector<double> cur_d_info = { d_cur, dd_cur, ddd_cur };
				vector<double> tgt_d_info = { target_d, target_dd, target_ddd };

				// generate trajectory
				int num_traj = ps.num_pred - num_prev_path;
				vector<double> s_traj = Trajectory::generateTrajectory(cur_s_info, tgt_s_info, 
																		num_traj, ps.pred_dt, ps.move_dt);
				vector<double> d_traj = Trajectory::generateTrajectory(cur_d_info, tgt_d_info, 
																		num_traj, ps.pred_dt, ps.move_dt);

				// use remaining points from previous path
				for (int i = 0; i < num_prev_path; ++i)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}

				// fill predicted xy trajectory
				for (int i = 0; i < num_traj; ++i)
				{
					vector<double> xy = getXY(s_traj[i], d_traj[i], waypoints_s, waypoints_x, waypoints_y);
					next_x_vals.push_back(xy[0]);
					next_y_vals.push_back(xy[1]);
				}
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// send to simulator
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
