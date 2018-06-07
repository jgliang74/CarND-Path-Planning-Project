#ifndef VEHICLE_H
#define VEHICLE_H
#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <string>

using namespace std;
int SPEED_LIMIT = 50;      // Maximum speed limit miles/hour
int BUFFER_V = 2;          // target speed = SPEED_LIMIT - BUFFER_V
float STOP_COST = 0.85;

class Vehicle {
public:

  map<string, int> lane_direction = {{"PLCL", -1}, {"LCL", -1}, {"LCR", 1}, {"PLCR", 1}};

  struct collider
  {
    bool collision ; // is there a collision?
    int  time;       // time collision happens
  };

  int    L = 1;
  int    preferred_buffer = 50; // impacts "keep lane" behavior.
  int    lane;
  double s;
  float  v;
  float  a;
  float  target_speed;
  int    lanes_available;
  float  max_acceleration;
  int    goal_lane;
  int    goal_s;
  //double ref_x;
  //double ref_y;
  //double ref_yaw;

  //vector<double> next_x_vals;
  //vector<double> next_y_vals;

  string state;

  /**
  * Constructor
  */
  Vehicle();
  Vehicle(int lane, float s, float v, float a, string state="CS");

  /**
  * Destructor
  */
  virtual ~Vehicle();

  vector<Vehicle> choose_next_state(map<int, vector<Vehicle>> predictions);
  vector<string>  successor_states();
  vector<Vehicle> generate_trajectory(string state, map<int, vector<Vehicle>> predictions);
  vector<float> get_kinematics(map<int, vector<Vehicle>> predictions, int lane);
  vector<Vehicle> constant_speed_trajectory();
  vector<Vehicle> keep_lane_trajectory(map<int, vector<Vehicle>> predictions);
  vector<Vehicle> lane_change_trajectory(string state, map<int, vector<Vehicle>> predictions);
  vector<Vehicle> prep_lane_change_trajectory(string state, map<int, vector<Vehicle>> predictions);
  void increment(int dt);
  float position_at(int t);
  bool get_vehicle_behind(map<int, vector<Vehicle>> predictions, int lane, Vehicle & rVehicle);
  bool get_vehicle_ahead(map<int, vector<Vehicle>> predictions, int lane, Vehicle & rVehicle);
  vector<Vehicle> generate_predictions(int horizon=2);
  void realize_next_state(vector<Vehicle> trajectory);
  void configure(vector<int> road_data);

};

#endif