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
#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>
#include <functional>
#include <string>
#include <iterator>
#include "spline.h"
#include "vehicle.h"
#include "cost.h"

using namespace std;

const float LANE_CHANGE = pow(10, 2);
const float SPEED_COST_WEIGTH  = pow(10, 5);

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

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
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


  int lane_num = 1;
  double ref_vel = 0.0;
  double MAX_ACCEL = 0.224;  // Maximum acceleration 10m/S^2
  int num_lanes = 3;
  int goal_s = int (max_s);
  int goal_lane = 1;

  vector<int> ego_config = {SPEED_LIMIT - BUFFER_V, num_lanes, INT_MAX, goal_lane};
  Vehicle ego = Vehicle(lane_num, 0, 0, 0);
  ego.configure(ego_config);
  ego.max_acceleration = MAX_ACCEL;
  ego.state = "KL";


  h.onMessage([&ego, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

                map<int, Vehicle> other_vehicles;
                map<int ,vector<Vehicle> > predictions;
           
                int prev_size = previous_path_x.size();

                ego.s = car_s;

                if(prev_size > 0)
                {
                   ego.s = car_s = end_path_s;
                }

                //
                // find other vehicle on the road using sensor fusion data
                //
                for(int i = 0; i < sensor_fusion.size(); i++)
                {
                   float  d = sensor_fusion[i][6];
                   int    l = (int)d/4;
                   double check_car_s = sensor_fusion[i][5];
                   double vx = sensor_fusion[i][3];
                   double vy = sensor_fusion[i][4];
                   double check_speed = sqrt(vx*vx+vy*vy);

                   check_car_s += ((double)prev_size*0.02*check_speed);

	           Vehicle vehicle = Vehicle(l, check_car_s, check_speed*2.24, 0);
		   vehicle.state = "CS";
		   other_vehicles.insert(std::pair<int, Vehicle>(i, vehicle));
		}

                map<int, Vehicle>::iterator it = other_vehicles.begin();
                while(it != other_vehicles.end())
                {
                   int v_id = it->first;
                   vector<Vehicle> preds = it->second.generate_predictions();
                   predictions[v_id] = preds;
                   it++;
                }

cout << "previous_path: " << prev_size << endl;
cout << "Vehicles on lane: " << ego.lane << endl;
cout << "Vehicles states: " << ego.state << endl;
cout << "Vehicles vel: " << ego.v << endl;
cout << "Vehicles s: " << ego.s << endl;

//Vehicle temp_vehicle1;
//Vehicle temp_vehicle2;
//if(ego.get_vehicle_ahead(predictions, ego.lane, temp_vehicle1)){
//   cout << "Ahead Vehicles s: " << temp_vehicle1.s << endl;
//}
//if(ego.get_vehicle_behind(predictions, ego.lane, temp_vehicle2)){
//cout << "Behind Vehicles s: " << temp_vehicle2.s << endl;
//}

//cout << endl;

        	vector<Vehicle> trajectory = ego.choose_next_state(predictions);
        	
                //
                // Make the state transition
                //

                ego.realize_next_state(trajectory);
cout << endl;
          	json msgJson;

                vector<double> ptsx;
                vector<double> ptsy;

                double ref_x = car_x;
                double ref_y = car_y;
                double ref_yaw = deg2rad(car_yaw);

                if(prev_size < 2)
                {
                   double prev_car_x = car_x - cos(car_yaw);
                   double prev_car_y = car_y - sin(car_yaw);

                   ptsx.push_back(prev_car_x);
                   ptsx.push_back(car_x);

                   ptsy.push_back(prev_car_y);
                   ptsy.push_back(car_y);
                }
                else
                {
                   ref_x = previous_path_x[prev_size-1];
                   ref_y = previous_path_y[prev_size-1];
               
                   double ref_x_prev = previous_path_x[prev_size-2];
                   double ref_y_prev = previous_path_y[prev_size-2];
                   ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);

                   ptsx.push_back(ref_x_prev);
                   ptsx.push_back(ref_x);

                   ptsy.push_back(ref_y_prev);
                   ptsy.push_back(ref_y);
                }

                vector<double> next_wp0 = getXY(car_s+30,(2+4*ego.lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
                vector<double> next_wp1 = getXY(car_s+60,(2+4*ego.lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
                vector<double> next_wp2 = getXY(car_s+90,(2+4*ego.lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);

                ptsx.push_back(next_wp0[0]);
                ptsx.push_back(next_wp1[0]);
                ptsx.push_back(next_wp2[0]);

                ptsy.push_back(next_wp0[1]);
                ptsy.push_back(next_wp1[1]);
                ptsy.push_back(next_wp2[1]);

                for(int i = 0; i< ptsx.size(); i++)
                {
                    double shift_x = ptsx[i] - ref_x;
                    double shift_y = ptsy[i] - ref_y;

                    ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
                    ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
                }

                tk::spline s;
                s.set_points(ptsx, ptsy);

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

                for(int i=0; i < previous_path_x.size(); i++)
                {
                    next_x_vals.push_back(previous_path_x[i]);
                    next_y_vals.push_back(previous_path_y[i]);
                }

                double target_x = 30.0;
                double target_y = s(target_x);
                double target_dist = sqrt((target_x*target_x) + (target_y*target_y));

                double x_add_on = 0;

                // starter code

                for(int i = 1; i <= 50 - previous_path_x.size(); i++)
                {
                    double N = (target_dist/(0.02*ego.v/2.24));
                    double x_point = x_add_on+(target_x)/N;
                    double y_point = s(x_point);

                    x_add_on = x_point;

                    double x_ref = x_point;
                    double y_ref = y_point;
                   
                    // rotate back to normal after rotating it earlier
                    x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
                    y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

                    x_point += ref_x;
                    y_point += ref_y;

                    next_x_vals.push_back(x_point);
                    next_y_vals.push_back(y_point);
                }

               
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

/**
The following is implementation of Vehicle class
 **/

/**
 * Initializes Vehicle
 */

Vehicle::Vehicle(){}

Vehicle::Vehicle(int lane, float s, float v, float a, string state) {

    this->lane = lane;
    this->s = s;
    this->v = v;
    this->a = a;
    this->state = state;
    max_acceleration = -1;

}

Vehicle::~Vehicle() {}


vector<Vehicle> Vehicle::choose_next_state(map<int, vector<Vehicle>> predictions) {
    /*
    Here you can implement the transition_function code from the Behavior Planning Pseudocode
    classroom concept. Your goal will be to return the best (lowest cost) trajectory corresponding
    to the next state.

    INPUT: A predictions map. This is a map of vehicle id keys with predicted
        vehicle trajectories as values. Trajectories are a vector of Vehicle objects representing
        the vehicle at the current timestep and one timestep in the future.
    OUTPUT: The the best (lowest cost) trajectory corresponding to the next ego vehicle state.

    */
    vector<string> states = successor_states();
    float cost;
    vector<float> costs;
    vector<string> final_states;
    vector<vector<Vehicle>> final_trajectories;

    for (vector<string>::iterator it = states.begin(); it != states.end(); ++it) {
cout << "State " << *it << endl;
        vector<Vehicle> trajectory = generate_trajectory(*it, predictions);
        if (trajectory.size() != 0) {
            cost = calculate_cost(*this, predictions, trajectory);
            costs.push_back(cost);
            final_trajectories.push_back(trajectory);

cout << "State " << *it << endl;
cout << "Cost: " << cost << endl;
        }
    }
//cout << endl;
    vector<float>::iterator best_cost = min_element(begin(costs), end(costs));
    int best_idx = distance(begin(costs), best_cost);
    return final_trajectories[best_idx];
}

vector<string> Vehicle::successor_states() {
    /*
    Provides the possible next states given the current state for the FSM 
    discussed in the course, with the exception that lane changes happen 
    instantaneously, so LCL and LCR can only transition back to KL.
    */
    vector<string> states;
    states.push_back("KL");

    string state = this->state;
    if(state.compare("KL") == 0) 
    {
        states.push_back("LCR");
        states.push_back("LCL");
    } 

#if 0
    string state = this->state;
    if(state.compare("KL") == 0) {
        states.push_back("PLCL");
        states.push_back("PLCR");
    } else if (state.compare("PLCL") == 0) {
        if (lane != 0) {
            states.push_back("PLCL");
            states.push_back("LCL");
        }
    } else if (state.compare("PLCR") == 0) {
        if (lane != lanes_available -1) {
            states.push_back("PLCR");
            states.push_back("LCR");
        }
    }
#endif
    //If state is "LCL" or "LCR", then just return "KL"
    return states;
}

vector<Vehicle> Vehicle::generate_trajectory(string state, map<int, vector<Vehicle>> predictions) {
    /*
    Given a possible next state, generate the appropriate trajectory to realize the next state.
    */
    vector<Vehicle> trajectory;
    if (state.compare("CS") == 0) {
        trajectory = constant_speed_trajectory();
    } else if (state.compare("KL") == 0) {
        trajectory = keep_lane_trajectory(predictions);
    } else if ((state.compare("LCL") == 0 && this->lane != 0)|| (state.compare("LCR") == 0 && this->lane != 2)) {
        trajectory = lane_change_trajectory(state, predictions);
    } else if ((state.compare("PLCL") == 0 && this->lane != 0)|| (state.compare("PLCR") == 0 && this->lane != 2)) {
        trajectory = prep_lane_change_trajectory(state, predictions);
    }
    return trajectory;
}

vector<float> Vehicle::get_kinematics(map<int, vector<Vehicle>> predictions, int lane) {
    /* 
    Gets next timestep kinematics (position, velocity, acceleration) 
    for a given lane. Tries to choose the maximum velocity and acceleration, 
    given other vehicle positions and accel/velocity constraints.
    */
    float max_velocity_accel_limit = this->max_acceleration + this->v;
    float new_position;
    float new_velocity;
    float new_accel;
    Vehicle vehicle_ahead;
    Vehicle vehicle_behind;

    if (get_vehicle_ahead(predictions, lane, vehicle_ahead)) 
    {

cout << "ahead vehicle velocity " << endl;
cout << "s: " << vehicle_ahead.s << endl;
cout << "v: " << vehicle_ahead.v << endl;

        if (get_vehicle_behind(predictions, lane, vehicle_behind))
        {

cout << "behind vehicle velocity " << endl;
cout << "v: " << vehicle_behind.v << endl;

            // new_velocity = vehicle_ahead.v; //must travel at the speed of traffic, regardless of preferred buffer
        } 
        else 
        {
            // float max_velocity_in_front = (vehicle_ahead.s - this->s - this->preferred_buffer) + vehicle_ahead.v - 0.5 * (this->a);
            // new_velocity = min(min(max_velocity_in_front, max_velocity_accel_limit), this->target_speed);
        }

        new_velocity = max(this->v - this->max_acceleration, this->max_acceleration);
    } 
    else
    {
        new_velocity = min(max_velocity_accel_limit, this->target_speed);
    }

    new_accel = new_velocity - this->v; //Equation: (v_1 - v_0)/t = acceleration
    new_position = this->s + new_velocity + new_accel / 2.0;

//cout << "keep lane " << endl;
//cout << "s: " << new_position << endl;
//cout << "v: " << new_velocity << endl;
//cout << "a: " << new_accel << endl;

    return{new_position, new_velocity, new_accel};
    
}

vector<Vehicle> Vehicle::constant_speed_trajectory() {
    /*
    Generate a constant speed trajectory.
    */
    float next_pos = position_at(1);
    vector<Vehicle> trajectory = {Vehicle(this->lane, this->s, this->v, this->a, this->state), 
                                  Vehicle(this->lane, next_pos, this->v, 0, this->state)};
    return trajectory;
}

vector<Vehicle> Vehicle::keep_lane_trajectory(map<int, vector<Vehicle>> predictions) {
    /*
    Generate a keep lane trajectory.
    */
    vector<Vehicle> trajectory = {Vehicle(this->lane, this->s, this->v, this->a, state)};
    vector<float> kinematics = get_kinematics(predictions, this->lane);
    float new_s = kinematics[0];
    float new_v = kinematics[1];
    float new_a = kinematics[2];

cout << "lane " << this->lane << endl;
cout << "s: " << kinematics[0] << endl;
cout << "v: " << kinematics[1] << endl;
cout << "a: " << kinematics[2] << endl;

    trajectory.push_back(Vehicle(this->lane, new_s, new_v, new_a, "KL"));
    return trajectory;
}

vector<Vehicle> Vehicle::prep_lane_change_trajectory(string state, map<int, vector<Vehicle>> predictions) {
    /*
    Generate a trajectory preparing for a lane change.
    */
    float new_s;
    float new_v;
    float new_a;
    Vehicle vehicle_behind;
    int new_lane = this->lane + lane_direction[state];
    vector<Vehicle> trajectory = {Vehicle(this->lane, this->s, this->v, this->a, this->state)};
    vector<float> curr_lane_new_kinematics = get_kinematics(predictions, this->lane);

    if (get_vehicle_behind(predictions, this->lane, vehicle_behind)) 
    {
        //Keep speed of current lane so as not to collide with car behind.
        new_s = curr_lane_new_kinematics[0];
        new_v = curr_lane_new_kinematics[1];
        new_a = curr_lane_new_kinematics[2];
        
    } 
    else
    {
        vector<float> best_kinematics;
        vector<float> next_lane_new_kinematics = get_kinematics(predictions, new_lane);
        //Choose kinematics with lowest velocity.
        if (next_lane_new_kinematics[1] < curr_lane_new_kinematics[1]) {
            best_kinematics = next_lane_new_kinematics;
        } else {
            best_kinematics = curr_lane_new_kinematics;
        }
        new_s = best_kinematics[0];
        new_v = best_kinematics[1];
        new_a = best_kinematics[2];
    }

    trajectory.push_back(Vehicle(this->lane, new_s, new_v, new_a, state));
    return trajectory;
}

vector<Vehicle> Vehicle::lane_change_trajectory(string state, map<int, vector<Vehicle>> predictions) {
    /*
    Generate a lane change trajectory.
    */
    int new_lane = this->lane + lane_direction[state];
    vector<Vehicle> trajectory;
    Vehicle next_lane_vehicle;
    //Check if a lane change is possible (check if another vehicle occupies that spot).
    for (map<int, vector<Vehicle>>::iterator it = predictions.begin(); it != predictions.end(); ++it) {
        next_lane_vehicle = it->second[0];
        if(next_lane_vehicle.lane == new_lane)
        {
           cout << "New lane distance " << abs(this->s - next_lane_vehicle.s) << endl;
           if ( abs(this->s - next_lane_vehicle.s) < 5) 
           {
            //If lane change is not possible, return empty trajectory.
            return trajectory;
           }
        }
    }
    trajectory.push_back(Vehicle(this->lane, this->s, this->v, this->a, this->state));
    vector<float> kinematics = get_kinematics(predictions, new_lane);

cout << "New lane " << new_lane << endl;
cout << "s: " << kinematics[0] << endl;
cout << "v: " << kinematics[1] << endl;
cout << "a: " << kinematics[2] << endl;

    trajectory.push_back(Vehicle(new_lane, kinematics[0], kinematics[1], kinematics[2], state));
    return trajectory;
}

void Vehicle::increment(int dt = 1) {
	this->s = position_at(dt);
}

float Vehicle::position_at(int t) {
    return this->s + this->v*t + this->a*t*t/2.0;
}

bool Vehicle::get_vehicle_behind(map<int, vector<Vehicle>> predictions, int lane, Vehicle & rVehicle) {
    /*
    Returns a true if a vehicle is found behind the current vehicle, false otherwise. The passed reference
    rVehicle is updated if a vehicle is found.
    */
    int max_s = -1;
    bool found_vehicle = false;
    Vehicle temp_vehicle;
    for (map<int, vector<Vehicle>>::iterator it = predictions.begin(); it != predictions.end(); ++it) {
        temp_vehicle = it->second[0];
        if (temp_vehicle.lane == lane && temp_vehicle.s < this->s && temp_vehicle.s > max_s) {
            max_s = temp_vehicle.s;
            rVehicle = temp_vehicle;
            found_vehicle = true;
        }
    }

    if((s - max_s < 30) && (s - max_s > 0))
    {
       found_vehicle = true;
    }

    return found_vehicle;
}

bool Vehicle::get_vehicle_ahead(map<int, vector<Vehicle>> predictions, int lane, Vehicle & rVehicle) {
    /*
    Returns a true if a vehicle is found ahead of the current vehicle, false otherwise. The passed reference
    rVehicle is updated if a vehicle is found.
    */
    int min_s = this->goal_s;
    bool found_vehicle = false;
    Vehicle temp_vehicle;
    for (map<int, vector<Vehicle>>::iterator it = predictions.begin(); it != predictions.end(); ++it) {
        temp_vehicle = it->second[0];
        if (temp_vehicle.lane == lane && temp_vehicle.s > this->s && temp_vehicle.s < min_s) {
            min_s = temp_vehicle.s;
            rVehicle = temp_vehicle;
            found_vehicle = true;
        }
    }

    if(min_s - s < 30)
    {
       found_vehicle = true;
    }
    else
    {
       found_vehicle = false;
    }

    return found_vehicle;
}

vector<Vehicle> Vehicle::generate_predictions(int horizon) 
{
    /*
    Generates predictions for non-ego vehicles to be used
    in trajectory generation for the ego vehicle.
    */

    vector<Vehicle> predictions;
    for(int i = 0; i < horizon; i++) 
    {
      float next_s = position_at(i);
      float next_v = 0;
      if (i < horizon-1) 
      {
        next_v = position_at(i+1) - s;
      }
      predictions.push_back(Vehicle(this->lane, next_s, next_v, 0));
    }

    return predictions;
}

void Vehicle::realize_next_state(vector<Vehicle> trajectory) {
    /*
    Sets state and kinematics for ego vehicle using the last state of the trajectory.
    */
    Vehicle next_state = trajectory[1];
    this->state = next_state.state;
    this->lane = next_state.lane;
    this->s = next_state.s;
    this->v = next_state.v;
    this->a = next_state.a;
}

void Vehicle::configure(vector<int> road_data) {
    /*
    Called by simulator before simulation begins. Sets various
    parameters which will impact the ego vehicle. 
    */
    target_speed = road_data[0];
    lanes_available = road_data[1];
    goal_s = road_data[2];
    goal_lane = road_data[3];
    // max_acceleration = road_data[4];
}

float goal_distance_cost(const Vehicle & vehicle, const vector<Vehicle> & trajectory, const map<int, vector<Vehicle>> & predictions, map<string, float> & data) {
    /*
    Cost increases based on distance of intended lane (for planning a lane change) and final lane of trajectory.
    Cost of being out of goal lane also becomes larger as vehicle approaches goal distance.
    */
    float cost;
    float distance = data["distance_to_goal"];
    if (distance > 0) {
        cost = 1 - 2*exp(-(abs(2.0*vehicle.goal_lane - data["intended_lane"] - data["final_lane"]) / distance));
    } else {
        cost = 1;
    }
    return cost;
}

float lane_change_cost(const Vehicle & vehicle, const vector<Vehicle> & trajectory, const map<int, vector<Vehicle>> & predictions, map<string, float> & data) {
    
    
    float cost;

    if(data["init_lane"] == data["final_lane"])
    {
       cost = 0;
    }
    else
    {
       cost = 1;
    }

    return cost;
}

float speed_cost(const Vehicle & vehicle, const vector<Vehicle> & trajectory, const map<int, vector<Vehicle>> & predictions, map<string, float> & data) {
    /*
    Cost becomes higher for trajectories with intended lane and final lane that have traffic slower than vehicle's target speed. 
    */

    //float proposed_speed_intended = lane_speed(predictions, data["intended_lane"]);
    //if (proposed_speed_intended < 0) {
    //    proposed_speed_intended = vehicle.target_speed;
    //}

    //float proposed_speed_final = lane_speed(predictions, data["final_lane"]);
    //if (proposed_speed_final < 0) {
    //    proposed_speed_final = vehicle.target_speed;
    //}
    
    //float cost = (2.0*vehicle.target_speed - proposed_speed_intended - proposed_speed_final)/vehicle.target_speed;
    float cost = 0;
    float proposed_speed_final = data["proposed_speed_final"];

    if(proposed_speed_final > SPEED_LIMIT)
    {
       cost = 1;
    }
    else if (proposed_speed_final < vehicle.target_speed)
    {
       cost = STOP_COST * ((vehicle.target_speed - proposed_speed_final)/vehicle.target_speed);
    }
    else
    {
       cost = (proposed_speed_final - vehicle.target_speed)/BUFFER_V;
    }

    return cost;
}

float lane_speed(const map<int, vector<Vehicle>> & predictions, int lane) {
    /*
    All non ego vehicles in a lane have the same speed, so to get the speed limit for a lane,
    we can just find one vehicle in that lane.
    */
    for (map<int, vector<Vehicle>>::const_iterator it = predictions.begin(); it != predictions.end(); ++it) {
        int key = it->first;
        Vehicle vehicle = it->second[0];
        if (vehicle.lane == lane && key != -1) {
            return vehicle.v;
        }
    }
    //Found no vehicle in the lane
    return -1.0;
}

float calculate_cost(const Vehicle & vehicle, const map<int, vector<Vehicle>> & predictions, const vector<Vehicle> & trajectory) {
    /*
    Sum weighted cost functions to get total cost for trajectory.
    */
    map<string, float> trajectory_data = get_helper_data(vehicle, trajectory, predictions);
    float cost = 0.0;

    //Add additional cost functions here.
    vector< function<float(const Vehicle & , const vector<Vehicle> &, const map<int, vector<Vehicle>> &, map<string, float> &)>> cf_list = {lane_change_cost, speed_cost};
    vector<float> weight_list = {LANE_CHANGE, SPEED_COST_WEIGTH};
    
    for (int i = 0; i < cf_list.size(); i++) {
        float new_cost = weight_list[i]*cf_list[i](vehicle, trajectory, predictions, trajectory_data);
        cout << "new Cost value: " << new_cost << endl;
        cost += new_cost;
    }

    return cost;

}

map<string, float> get_helper_data(const Vehicle & vehicle, const vector<Vehicle> & trajectory, const map<int, vector<Vehicle>> & predictions) {
    /*
    Generate helper data to use in cost functions:
    indended_lane: the current lane +/- 1 if vehicle is planning or executing a lane change.
    final_lane: the lane of the vehicle at the end of the trajectory.
    distance_to_goal: the distance of the vehicle to the goal.

    Note that indended_lane and final_lane are both included to help differentiate between planning and executing
    a lane change in the cost functions.
    */
    map<string, float> trajectory_data;
    Vehicle trajectory_first = trajectory[0];
    Vehicle trajectory_last = trajectory[1];
    float intended_lane;

    if (trajectory_last.state.compare("PLCL") == 0) {
        intended_lane = trajectory_last.lane + 1;
    } else if (trajectory_last.state.compare("PLCR") == 0) {
        intended_lane = trajectory_last.lane - 1;
    } else {
        intended_lane = trajectory_last.lane;
    }

    float distance_to_goal = vehicle.goal_s - trajectory_last.s;
    float final_lane = trajectory_last.lane;
    trajectory_data["intended_lane"] = intended_lane;
    trajectory_data["final_lane"] = final_lane;
    trajectory_data["distance_to_goal"] = distance_to_goal;
    trajectory_data["init_lane"] = trajectory_first.lane;
    trajectory_data["proposed_speed_final"] = trajectory_last.v;
    return trajectory_data;
}