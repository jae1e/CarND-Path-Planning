#include <fstream>
#define _USE_MATH_DEFINES
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
double mph2mps(double x) { return x * 0.44704; }

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
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y, const vector<double> &maps_s)
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
	double frenet_s = maps_s[0];
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenetClosest(double x, double y, int closest_wp,
								const vector<double> &maps_x, const vector<double> &maps_y, const vector<double> &maps_s)
{
	int prev_wp = closest_wp-1;
	if(closest_wp == 0)
	{
		throw;
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[closest_wp]-maps_x[prev_wp];
	double n_y = maps_y[closest_wp]-maps_y[prev_wp];
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
	double frenet_s = maps_s[0];
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

enum BehaviorStatus
{
	Start = -1,
	KeepLane = 0,
	ChangeLane = 1
};

enum LaneStatus
{
	Busy = 0,
	Free = 1
};

struct ParameterSet
{
public: // general parameters
	double len_track_s = 6945.554;

	double max_s = 0;
	
	int min_traj = 20;
	
	double iter_dt = 0.02;
		
	double max_accel = 8.0;

	double max_jerk = 5.0;
	
	double max_speed = 22.2; // 50 MPH = 22.352 m/s;
	
	double lane_width = 4.0;
	
	int num_lanes = 3;
	
	double accel_coeff = 0.5; // velocity quickly converges when the value goes high
	
public: // start parameters
	double start_duration = 5.0;
	
	double start_speed = 12.0;
	
public: // lane keep parameters
	double lane_keep_duration = 0.5;
	
	double lane_keep_speed_change = 0.025;

	double max_d_comp_curve_angle = pi() / 8; // maximum compensation curve angle

	double lane_curve_duration = 1.0;

	double ds_curve_comp_coeff = 5.0;

	double d_curve_comp_coeff = 0.8;

	double max_d_deviation = 0.9;
		
public: // lane change parameters
	double lane_change_duration = 3.0;
	
	double lane_change_distance = 50.0;

	double lane_change_speed_decay = 0.05;
	
	double lane_change_cost_coeff = 1.1;

	double lane_change_front_safety_distance = 25.0;

	double lane_change_back_safety_distance = 20.0;
	
public: // safety parameters
	double safety_change_duration = 2.0;

	double safety_distance = 20.0;

	double safety_speed_change = 0.2;

public: // interpoaltion parameters
	int num_src_waypoints = 10;
	
	double interpolation_interval = 0.5;
	
public: // status parameters
	int num_cycle = 0;

	double prev_car_s_mod = -1;

	int current_lane_id = 1;
	
	int target_lane_id = 1;
	
	double last_traj_s = 0;
	
	double last_traj_ds = 0;
	
	double last_traj_dds = 0;
	
	double last_traj_d = 0;
	
	double last_traj_dd = 0;
	
	double last_traj_ddd = 0;
	
	BehaviorStatus current_status = BehaviorStatus::Start;
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

	bool crazymode = false;
	if (crazymode)
	{
		ps.max_speed = 120;
		ps.start_speed = 20.0;
		ps.lane_keep_speed_change = 0.3;
		ps.lane_change_front_safety_distance = 1000000000000;
		ps.lane_change_back_safety_distance = 1000000000000;
		ps.lane_change_distance = 0;
		ps.safety_distance = 0;
		ps.num_src_waypoints = 20;
	}

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

          	double car_s_mod = j[1]["s"];
			if (car_s_mod < ps.prev_car_s_mod - 0.5 * ps.len_track_s)
			{
				ps.num_cycle += 1;
			}
			ps.prev_car_s_mod = car_s_mod;
			double car_s = car_s_mod + ps.num_cycle * ps.len_track_s;

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

			printf("-------------------------------------------\n");
			printf("\n");
			printf("------------------\n");
			printf(" current car info\n");
			printf("------------------\n");
			printf("s: %.2f\t d: %.2f\n", car_s, car_d);
			printf("x: %.2f\t y: %.2f\n", car_x, car_y);
			printf("speed: %.2f\t yaw: %.2f\n", mph2mps(car_speed), deg2rad(car_yaw));

          	//////////////////////////////////////////////////////////////////////////////////////////////////////
			// interpolate waypoints
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			vector<double> src_waypoints_s, src_waypoints_x, src_waypoints_y, src_waypoints_dx, src_waypoints_dy;
			vector<double> waypoints_s, waypoints_x, waypoints_y, waypoints_dx, waypoints_dy;
			
			{
				int num_total_waypoints = map_waypoints_x.size();
				int next_waypoint_id = NextWaypoint(car_x, car_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);
				
				int begin_waypoint_id = next_waypoint_id - 0.5 * ps.num_src_waypoints;
				int end_waypoint_id = next_waypoint_id + 0.5 *  ps.num_src_waypoints;
				
				for (int wid = begin_waypoint_id; wid < end_waypoint_id; ++wid)
				{
					int wid_mod = wid;
					while (wid_mod > num_total_waypoints)
					{
						wid_mod -= num_total_waypoints;
					}
					while (wid_mod < 0)
					{
						wid_mod += num_total_waypoints;
					}

					double sval = map_waypoints_s[wid_mod];
					while (sval > car_s + 0.5 * ps.len_track_s)
					{
						sval -= ps.len_track_s;
					}
					while (sval < car_s - 0.5 * ps.len_track_s)
					{
						sval += ps.len_track_s;
					}

					if (src_waypoints_s.empty() || sval > src_waypoints_s.back() + 1e-2)
					{
						src_waypoints_s.push_back(sval);

						src_waypoints_x.push_back(map_waypoints_x[wid_mod]);
						src_waypoints_y.push_back(map_waypoints_y[wid_mod]);
						src_waypoints_dx.push_back(map_waypoints_dx[wid_mod]);
						src_waypoints_dy.push_back(map_waypoints_dy[wid_mod]);
					}
				}

				for (int i = 0; i < src_waypoints_s.size() - 1; ++i)
				{
					double p1 = src_waypoints_s[i];
					double p2 = src_waypoints_s[i+1];
					int num_points = (int)ceil((p2 - p1) / ps.interpolation_interval);
		
					for (int j = 0; j < num_points; ++j)
					{
						waypoints_s.push_back(p1 + j * ps.interpolation_interval);
					}
				}
				waypoints_s.push_back(src_waypoints_s.back());

				//for (int i = 0; i != src_waypoints_s.size(); ++i)
				//{
				//	printf("s: %f x: %f y: %f\n", src_waypoints_s[i], src_waypoints_x[i], src_waypoints_y[i]);
				//}

				Utils::interpolatePoints(src_waypoints_s, src_waypoints_x, ps.interpolation_interval, waypoints_x);
				Utils::interpolatePoints(src_waypoints_s, src_waypoints_y, ps.interpolation_interval, waypoints_y);
				Utils::interpolatePoints(src_waypoints_s, src_waypoints_dx, ps.interpolation_interval, waypoints_dx);
				Utils::interpolatePoints(src_waypoints_s, src_waypoints_dy, ps.interpolation_interval, waypoints_dy);

				ps.max_s = waypoints_s.back() + ps.len_track_s;
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// check lane status
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			double car_ds = 0;
			
			vector<double> lane_preceding_s_dist(ps.num_lanes, ps.max_s);
			vector<double> lane_preceding_ds(ps.num_lanes, 0);
			vector<LaneStatus> lane_status(ps.num_lanes, LaneStatus::Free);
			
			{
				vector<double> lane_pred_s(ps.num_lanes, ps.max_s);
				
				{
					int closest_wp = ClosestWaypoint(car_x, car_y, waypoints_x, waypoints_y);
					double waypoint_angle = atan2(waypoints_y[closest_wp + 1] - waypoints_y[closest_wp],
												  waypoints_x[closest_wp + 1] - waypoints_x[closest_wp]);
					double frenet_angle = Utils::frenetAngle(deg2rad(car_yaw), waypoint_angle);
					car_ds = mph2mps(car_speed) * cos(frenet_angle);
				}
				
				// predicted s of lane after lane change time
				vector<double> lane_pred_min_s(ps.num_lanes, ps.max_s);
				double pred_time = 0.5 * ps.lane_keep_duration + ps.lane_change_duration;
				double pred_ego_s = car_s + car_ds * pred_time;
				double pred_ego_lc_s = car_s + car_ds * (1 - ps.lane_change_speed_decay) * pred_time;
				
				// fill ds and dspeed of each lane
				for (auto sfit = sensor_fusion.begin(); sfit != sensor_fusion.end(); ++sfit)
				{
					// int id = sfit->at(0);
					double x = sfit->at(1), y = sfit->at(2);
					double vx = sfit->at(3), vy = sfit->at(4);
					double v = sqrt(vx * vx + vy * vy);
					double s = (double)sfit->at(5) + ps.num_cycle * ps.len_track_s;
					double d = sfit->at(6);
					double yaw = atan2(vy, vx);
					
					int closest_wp = ClosestWaypoint(x, y, waypoints_x, waypoints_y);
					double waypoint_angle = atan2(waypoints_y[closest_wp + 1] - waypoints_y[closest_wp],
												  waypoints_x[closest_wp + 1] - waypoints_x[closest_wp]);
					double frenet_angle = Utils::frenetAngle(yaw, waypoint_angle);
					double ds = v * cos(frenet_angle);
					double dd = v * sin(frenet_angle);
					// printf("other car's v ds dd: %f %f %f\n", v, ds, dd);
					// printf("other car's yaw wpa frenet: %f %f %f\n", yaw, waypoint_angle, frenet_angle);
					
					// current info
					int lane_id = (int)(d / ps.lane_width);
					if (lane_id > -1 && lane_id < ps.num_lanes
						&& (s > car_s && s - car_s < lane_preceding_s_dist[lane_id]))
					{
						lane_preceding_s_dist[lane_id] = s - car_s;
						lane_preceding_ds[lane_id] = ds;
					}
					
					// check lane busy with current info and predicted info
					if (lane_id > -1 && lane_id < ps.num_lanes
						&& (s < car_s + ps.lane_change_front_safety_distance
							&& s > car_s - ps.lane_change_back_safety_distance))
					{
						lane_status[lane_id] = LaneStatus::Busy;
					}
					
					double pred_s = s + ds * pred_time;
					double pred_d = d + dd * pred_time;
					int pred_lane_id = (int)(pred_d / ps.lane_width);
					if (pred_lane_id > -1 && pred_lane_id < ps.num_lanes
						&& ((pred_s < pred_ego_s + ps.lane_change_front_safety_distance
							&& pred_s > pred_ego_s - ps.lane_change_back_safety_distance)
							|| (pred_s < pred_ego_lc_s + ps.lane_change_front_safety_distance
								&& pred_s > pred_ego_lc_s - ps.lane_change_back_safety_distance)))
					{
						lane_status[pred_lane_id] = LaneStatus::Busy;
					}
				}
			}
			
			vector<string> lss(ps.num_lanes);
			for (int i = 0; i < ps.num_lanes; ++i)
			{
				lss[i] = lane_status[i] == LaneStatus::Busy ? "busy" : "free";
			}
			
			printf("\n");
			printf("-----------\n");
			printf(" lane info\n");
			printf("-----------\n");
			printf("0: %s\t 1: %s\t 2: %s\n", lss[0].c_str(), lss[1].c_str(), lss[2].c_str());
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// make decision
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
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
					// if current lane's min distance and delta speed is not enough,
					if (lane_preceding_s_dist[ps.current_lane_id] < ps.lane_change_distance + ps.lane_change_front_safety_distance
						&& lane_preceding_ds[ps.current_lane_id] < car_ds)
					{
						// calculate lane cost
						double current_pred = Utils::prececedingDistPred(lane_preceding_s_dist[ps.current_lane_id],
																		 lane_preceding_ds[ps.current_lane_id],
																		 ps.lane_keep_duration);
						
						double left_pred = 0, right_pred = 0;
						if (ps.current_lane_id > 0
							&& lane_status[ps.current_lane_id - 1] == LaneStatus::Free)
						{
							left_pred = Utils::prececedingDistPred(lane_preceding_s_dist[ps.current_lane_id - 1],
																   lane_preceding_ds[ps.current_lane_id - 1],
																   ps.lane_keep_duration);
						}
						if (ps.current_lane_id < ps.num_lanes - 1
							&& lane_status[ps.current_lane_id + 1] == LaneStatus::Free)
						{
							right_pred = Utils::prececedingDistPred(lane_preceding_s_dist[ps.current_lane_id + 1],
																	lane_preceding_ds[ps.current_lane_id + 1],
																	ps.lane_keep_duration);
						}
						
						// make decision
						double max_pred = max(left_pred, right_pred);
						if (max_pred > current_pred * ps.lane_change_cost_coeff)
						{
							ps.target_lane_id = left_pred < right_pred ? ps.current_lane_id + 1 : ps.current_lane_id - 1;
							ps.current_status = BehaviorStatus::ChangeLane;
						}
					}
				}
			}
			
			std::string decision = ps.current_status == BehaviorStatus::ChangeLane ? "change lane" : "keep lane";
			
			printf("\n");
			printf("--------------------\n");
			printf(" lane decision info\n");
			printf("--------------------\n");
			printf("decision: %s\n", decision.c_str());
			printf("current lane id: %d\n", ps.current_lane_id);
			printf("target lane id: %d\n", ps.target_lane_id);
			printf("target lane s dist: %f\n", lane_preceding_s_dist[ps.target_lane_id]);
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// target generation
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			
			double cur_s = ps.last_traj_s == 0 ? car_s : ps.last_traj_s;
			double cur_ds = ps.last_traj_ds;
			double cur_dds = ps.last_traj_dds;
			double cur_d = ps.last_traj_s == 0 ? car_d : ps.last_traj_d;
			double cur_dd = ps.last_traj_dd;
			double cur_ddd = ps.last_traj_ddd;
			
			double target_s = 0, target_ds = 0, target_dds = 0;
			double target_d = 0, target_dd = 0, target_ddd = 0;
			
			double pred_duration = 0;
			
			{
				// // conversion test
				// if (previous_path_x.size() > ps.num_emergency_traj_keep)
				// {
				// 	double x = previous_path_x[ps.num_emergency_traj_keep - 1];
				// 	double y = previous_path_y[ps.num_emergency_traj_keep - 1];
					
				// 	vector<double> sd = getFrenet(x, y, deg2rad(car_yaw), waypoints_x, waypoints_y, waypoints_s);
					
				// 	int closest_wp = ClosestWaypoint(x, y, waypoints_x, waypoints_y);
				// 	double waypoint_angle = atan2(waypoints_y[closest_wp + 1] - waypoints_y[closest_wp],
				// 								  waypoints_x[closest_wp + 1] - waypoints_x[closest_wp]);
				// 	double frenet_angle = deg2rad(car_yaw) - waypoint_angle;
				// 	double ds = car_speed * cos(frenet_angle);
				// 	double dd = car_speed * sin(frenet_angle);

				// 	printf("conversion x y angle s d: %f %f %f %f %f\n", x, y, waypoint_angle, sd[0], sd[1]);
				// }

				// calculate target d info
				target_ddd = 0;
				target_dd = 0;
				target_d = (0.5 + ps.target_lane_id) * ps.lane_width;

				bool ds_decreased = false;

				// calculate target s info
				if (ps.current_status == BehaviorStatus::Start)
				{
					pred_duration = ps.start_duration;
					target_dds = 0;
					target_ds = ps.start_speed;

					printf("status: starting\n");
				}
				else if (ps.current_status == BehaviorStatus::ChangeLane)
				{
					pred_duration = ps.lane_change_duration;
					target_dds = 0;
					target_ds = cur_ds * (1.0 - ps.lane_change_speed_decay);
					ds_decreased = true;

					printf("status: changing lane\n");
				}
				else // if (ps.current_status == BehaviorStatus::KeepLane)
				{
					target_dds = 0;
					
					// emergency
					if (lane_preceding_s_dist[ps.target_lane_id] < ps.safety_distance)
					{
						pred_duration = ps.safety_change_duration;

						double prec_ds = lane_preceding_ds[ps.target_lane_id];
						target_ds = (1 - ps.safety_speed_change) * cur_ds;
						ds_decreased = true;

						// printf("target_ds: %f\n", target_ds);
						printf("status: car detected within safety range\n");
					}
					// match speed to preceding car
					else if (lane_preceding_s_dist[ps.target_lane_id] < ps.lane_change_distance)
					{
						pred_duration = ps.lane_keep_duration;

						double prec_ds = lane_preceding_ds[ps.target_lane_id];
						if (prec_ds > cur_ds)
						{
							target_ds = min(ps.max_speed, 
											min(prec_ds, (1 + ps.lane_keep_speed_change) * cur_ds));
						}
						else
						{
							target_ds = max(prec_ds, (1 - ps.lane_keep_speed_change) * cur_ds);
							ds_decreased = true;
						}

						ds_decreased = true;

						printf("status: need lane change\n");
					}
					// match speed to maximum speed
					else
					{
						pred_duration = ps.lane_keep_duration;

						double target_speed = ps.max_speed;
						if (target_speed > cur_ds)
						{
							target_ds = min(ps.max_speed, 
											min(target_speed, (1 + ps.lane_keep_speed_change) * cur_ds));
						}
						else
						{
							target_ds = max(target_speed, (1 - ps.lane_keep_speed_change) * cur_ds);
							ds_decreased = true;
						}

						printf("status: free to go\n");
					}
				}

				// curve compensation
				if (previous_path_x.size() > 1)
				{
					double x0 = previous_path_x.front();
					double y0 = previous_path_y.front();
					double x1 = previous_path_x.back();
					double y1 = previous_path_y.back();
					int wp0 = ClosestWaypoint(x0, y0, waypoints_x, waypoints_y);
					int wp1 = ClosestWaypoint(x1, y1, waypoints_x, waypoints_y);
					double angle0 = atan2(waypoints_y[wp0 + 1] - waypoints_y[wp0],
											waypoints_x[wp0 + 1] - waypoints_x[wp0]);
					double angle1 = atan2(waypoints_y[wp1 + 1] - waypoints_y[wp1],
											waypoints_x[wp1 + 1] - waypoints_x[wp1]);
					double angle_diff = Utils::normalizeAngle(angle0 - angle1);
					// printf("angle diff: %f\n", angle_diff);

					double target_lane_d = (0.5 + ps.current_lane_id) * ps.lane_width;
					// prevent outside lane
					if (ps.current_status == BehaviorStatus::KeepLane 
						&& ((ps.target_lane_id == 2 && cur_d - target_lane_d > ps.max_d_deviation)
							|| (ps.target_lane_id == 0 && target_lane_d - cur_d > ps.max_d_deviation)))
					{
						pred_duration = max(pred_duration, ps.lane_curve_duration);

						target_d = target_lane_d + 0.5 * (target_lane_d - cur_d);
					}
					// target d compensation according to angle difference, not to go out of the lane
					else if (abs(cur_d - target_lane_d) > 0.5 * ps.max_d_deviation)
					{
						pred_duration = max(pred_duration, ps.lane_curve_duration);

						double angle_sign = angle_diff > 0 ? 1 : -1;
						double d_comp_angle = abs(angle_diff) < ps.max_d_comp_curve_angle ? abs(angle_diff) : ps.max_d_comp_curve_angle;
						d_comp_angle *= angle_sign;
						target_d += ps.d_curve_comp_coeff * (d_comp_angle / ps.max_d_comp_curve_angle) * ps.lane_width;
					}

					if (!ds_decreased)
					{
						double ds_comp = cos(ps.ds_curve_comp_coeff * angle_diff);
						//printf("target ds comp: %f %f \n", target_ds, ds_comp);
						target_ds = target_ds * ds_comp;
					}
				}

				// calculate efficient target s
				target_s = Utils::efficientDistance(cur_s, cur_ds, target_ds, pred_duration, ps.accel_coeff);
			}
			
//			printf("----- cur s, ds, dds: %f %f %f\n", cur_s, cur_ds, cur_dds);
//			printf("----- target s, ds, dds: %f %f %f\n", target_s, target_ds, target_dds);
//			printf("----- cur d, dd, ddd: %f %f %f\n", cur_d, cur_dd, cur_ddd);
//			printf("----- target d, dd, ddd: %f %f %f\n", target_d, target_dd, target_ddd);
//			printf("----- pred_duration: %f\n", pred_duration);
			
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// generate trajectory
			//////////////////////////////////////////////////////////////////////////////////////////////////////
						
			{
				int num_prev_path = previous_path_x.size();
				int num_traj_calc = pred_duration / ps.iter_dt;
				int num_traj_fill = num_traj_calc;// ps.current_status == BehaviorStatus::Start ? num_traj_calc : ps.min_traj - num_prev_path;
				// int num_traj_fill = num_traj_calc;

				// use remaining points from previous path
				for (int i = 0; i < num_prev_path; ++i)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
					
					// printf("old traj xy: %f %f\n", (double)previous_path_x[i], (double)previous_path_y[i]);
				}
				
				if (num_prev_path < ps.min_traj)
				{
					vector<double> cur_s_info = { cur_s, cur_ds, cur_dds };
					vector<double> tgt_s_info = { target_s, target_ds, target_dds };
					vector<double> cur_d_info = { cur_d, cur_dd, cur_ddd };
					vector<double> tgt_d_info = { target_d, target_dd, target_ddd };

					printf("cur_s_info\n");
					for (int i = 0; i < 3; i++)
					{
						printf("%f\n", cur_s_info[i]);
					}
					printf("tgt_s_info\n");
					for (int i = 0; i < 3; i++)
					{
						printf("%f\n", tgt_s_info[i]);
					}
					
					vector<double> traj_s = Utils::jerkMinimizingTrajectory(cur_s_info, tgt_s_info,
																			num_traj_calc, pred_duration, ps.iter_dt);
					vector<double> traj_d = Utils::jerkMinimizingTrajectory(cur_d_info, tgt_d_info,
																			num_traj_calc, pred_duration, ps.iter_dt);
					
					//for (int i = 0; i < num_prev_path; ++i)
					//{
					//	double x = previous_path_x[i];
					//	double y = previous_path_y[i];
					//	int closest_wp = ClosestWaypoint(x, y, waypoints_x, waypoints_y);
					//	double waypoint_angle = atan2(waypoints_y[closest_wp + 1] - waypoints_y[closest_wp],
					//								  waypoints_x[closest_wp + 1] - waypoints_x[closest_wp]);
					//	
					//	vector<double> sd = getFrenet(x, y, waypoint_angle, waypoints_x, waypoints_y, waypoints_s);
					//	
					//	printf("old traj sd: %f %f\n", sd[0], sd[1]);
					//	printf("old traj xy: %f %f\n", (double)x, (double)y);
					//}

					// add new points to trajectory
					double pt_x = next_x_vals.empty() ? car_x : next_x_vals.back();
					double pt_y = next_y_vals.empty() ? car_y : next_y_vals.back();

					double max_dist = ps.max_speed * ps.iter_dt;

					for (int i = 1; i < num_traj_fill; ++i)
					{
						double del_s = traj_s[i] - traj_s[i - 1];
						double del_d = traj_d[i] - traj_d[i - 1];

						// prevent exceed speed limit
						double del_dist = sqrt(del_s * del_s + del_d * del_d);
						if (del_dist > max_dist)
						{
							del_s *= max_dist / del_dist;
							del_d *= max_dist / del_dist;
						}

						int last_wp = ClosestWaypoint(pt_x, pt_y, waypoints_x, waypoints_y);

						double dx = waypoints_dx[last_wp - 1];
						double dy = waypoints_dy[last_wp - 1];
						double sx = -dy;
						double sy = dx;

						pt_x += del_s * sx + del_d * dx;
						pt_y += del_s * sy + del_d * dy;

						next_x_vals.push_back(pt_x);
						next_y_vals.push_back(pt_y);

						// printf("new traj sd: %f %f\n", traj_s[i], traj_d[i]);
						// printf("new traj xy: %f %f\n", xy[0], xy[1]);
					}
				}

				int last_id = next_x_vals.size() - 1;

				double last_x3 = next_x_vals[last_id];
				double last_y3 = next_y_vals[last_id];
				double last_x2 = next_x_vals[last_id - 1];
				double last_y2 = next_y_vals[last_id - 1];
				double last_x1 = next_x_vals[last_id - 2];
				double last_y1 = next_y_vals[last_id - 2];

				double angle = atan2(last_y3 - last_y2, last_x3 - last_x2);

				// int last_wp = NextWaypoint(last_x3, last_y3, angle, waypoints_x, waypoints_y);
				int last_wp = ClosestWaypoint(last_x3, last_y3, waypoints_x, waypoints_y);

				vector<double> last_sd = getFrenetClosest(last_x3, last_y3, last_wp,
															waypoints_x, waypoints_y, waypoints_s);

				double dx = waypoints_dx[last_wp - 1];
				double dy = waypoints_dy[last_wp - 1];
				double sx = -dy;
				double sy = dx;

				double vx2 = (last_x3 - last_x2) / ps.iter_dt;
				double vy2 = (last_y3 - last_y2) / ps.iter_dt;
				double vx1 = (last_x2 - last_x1) / ps.iter_dt;
				double vy1 = (last_y2 - last_y1) / ps.iter_dt;

				double ax = (vx2 - vx1) / ps.iter_dt;
				double ay = (vy2 - vy1) / ps.iter_dt;

				// save trajectory last point
				ps.last_traj_s = last_sd[0];
				ps.last_traj_ds = vx2 * sx + vy2 * sy;
				ps.last_traj_dds = ax * sx + ay * sy;
				ps.last_traj_d = last_sd[1];
				ps.last_traj_dd = vx2 * dx + vy2 * dy;
				ps.last_traj_ddd = ax * dx + ay * dy;
				
				printf("\n");
				printf("-----------------\n");
				printf(" trajectory info\n");
				printf("-----------------\n");
				printf("target s: %.2f\t d: %.2f\n", ps.last_traj_s, ps.last_traj_d);
				printf("target ds: %.2f\t dd: %.2f\n", ps.last_traj_ds, ps.last_traj_dd);
				printf("target dds: %.2f\t ddd: %.2f\n", ps.last_traj_dds, ps.last_traj_ddd);
				printf("remaining count: %d\n", num_prev_path);
				printf("calculation count: %d\n", num_traj_calc);
				printf("new count: %d\n", num_traj_fill);
				printf("\n");
				printf("-------------------------------------------\n");
			}
			
			{
				// change status from start to keep lane
				if (ps.current_status == BehaviorStatus::Start)
				{
					ps.current_status = BehaviorStatus::KeepLane;
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
    // std::cout << "Connected!!!" << std::endl;
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
