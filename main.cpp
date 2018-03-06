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
#include "spline.h"

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
  int counter;
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
  
  double slow_car = 0.0;
  int lane = 0; //middle lane walkthru video, was 1
  double ref_vel = 0.0; //max mph speed walkthru video
  double check_speed = 0.0;
  
  double cost_ahead_left = 0; //was 3.0
  double cost_behind_left = 0;
  double cost_ahead_middle = 0;
  double cost_behind_middle = 0;
  double cost_ahead_right = 0;
  double cost_behind_right = 0;
  
  double cost_benefit_left_ahead = 0; //was 1000
  double cost_benefit_left_behind = 0;
  double cost_benefit_right_ahead = 0;
  double cost_benefit_right_behind = 0;
  double cost_benefit_middle_ahead = 0;
  double cost_benefit_middle_behind = 0;
  
  double no_cars_left_ahead = 0;
  double no_cars_left_behind = 0;
  double no_cars_middle_ahead = 0;
  double no_cars_middle_behind = 0;
  double no_cars_right_ahead = 0;
  double no_cars_right_behind = 0;
  
  double car_speed_sum = 0.0;
  double ave_car_speed = 0.0;
  int k = 1;
  int m = 0;
  int cycle = 0;
  h.onMessage([&cycle,&m,&cost_ahead_middle,&cost_behind_middle,&cost_ahead_right,&cost_behind_right,&cost_ahead_left,&cost_behind_left,&no_cars_left_ahead,&no_cars_left_behind,&no_cars_middle_ahead,&no_cars_middle_behind,&no_cars_right_ahead,&no_cars_right_behind,&k,&car_speed_sum,&ave_car_speed,&slow_car,&cost_benefit_left_ahead,&cost_benefit_middle_ahead,&cost_benefit_middle_behind,&cost_benefit_left_behind,&cost_benefit_right_ahead,&cost_benefit_right_behind,&check_speed,&ref_vel,&lane,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
	{

      auto s = hasData(data);

		if (s != "")
		{
			auto j = json::parse(s);
        
			string event = j[0].get<string>();
        
			if (event == "telemetry")
			{
				// j[1] is the data JSON object
          
				// Main car's localization Data
				double car_x = j[1]["x"];
				double car_y = j[1]["y"];
				double car_s = j[1]["s"];
				double car_d = j[1]["d"];
				double car_yaw = j[1]["yaw"];
				double car_speed = j[1]["speed"];
				
				
				
				if(k<=40)
				{
					car_speed_sum = car_speed_sum + car_speed;
					ave_car_speed = car_speed_sum/k;
					k = k + 1;
				}
				else
				{
					car_speed_sum = 0; //reset
					k = 1; //reset
				}
				
				//determine the car lane
				if(car_d > 0 & car_d < 4) //4 meters is width of each lane
				{
					lane = 0; //left lane
				}
				else if (car_d > 4 & car_d < 8)
				{
					lane = 1; //middle lane
				}
				else if (car_d > 8 & car_d < 12)
				{
					lane = 2; //right lane
				}
				//cout<<"car lane"<<lane<<endl; //lane id for car working correctly
			
				// Previous path data given to the Planner
				auto previous_path_x = j[1]["previous_path_x"];
				auto previous_path_y = j[1]["previous_path_y"];
				// Previous path's end s and d values 
				double end_path_s = j[1]["end_path_s"];
				double end_path_d = j[1]["end_path_d"];

				// Sensor Fusion Data, a list of all other cars on the same side of the road.
				auto sensor_fusion = j[1]["sensor_fusion"];

				int prev_size = previous_path_x.size();
			
				if(prev_size > 0)
				{
					car_s = end_path_s;
				}
				
				bool too_close = false;
				
				if(cycle < 30) //initialization for 30 cycles, hold car from lane changes for sensor updates
				{
					no_cars_left_ahead = 0; 
					no_cars_left_behind = 0; 
					no_cars_middle_ahead = 0; 
					no_cars_middle_behind = 0; 
					no_cars_right_ahead = 0; 
					no_cars_right_behind = 0; 
				
					cost_ahead_left = 0; 
					cost_behind_left = 0;
					cost_ahead_middle = 0;
					cost_behind_middle = 0;
					cost_ahead_right = 0;
					cost_behind_right = 0;
				
					//minimum costs for nearest cars
					cost_benefit_left_ahead = 0; //middle lane usage
					cost_benefit_right_ahead = 0; //middle lane usage
					cost_benefit_left_behind = 0; //middle lane usage
					cost_benefit_right_behind = 0; //middle lane usage
					cost_benefit_middle_ahead = 0; //left and right lane usage
					cost_benefit_middle_behind = 0; //left and right lane usage
				}
			
				else
				{
					for(m = 0; m < 1; m++) //initialization of costs after each sensor sweep
					{
						no_cars_left_ahead = 0; //reset
						no_cars_left_behind = 0; //reset
						no_cars_middle_ahead = 0; //reset
						no_cars_middle_behind = 0; //reset
						no_cars_right_ahead = 0; //reset
						no_cars_right_behind = 0; //reset
				
						//minimum costs for nearest cars
						cost_benefit_left_ahead = 7000.0; //reset for middle lane usage
						cost_benefit_right_ahead = 7000.0; //reset for middle lane usage
						cost_benefit_left_behind = 7000.0; //reset for middle lane usage
						cost_benefit_right_behind = 7000.0; //reset for middle lane usage
						cost_benefit_middle_ahead = 7000.0; //reset for left and right lane usage
						cost_benefit_middle_behind = 7000.0; //reset for left and right lane usage
					}
				}
				
				cycle = cycle + 1; //advance counter for cycle
				
				for(int i = 0; i < sensor_fusion.size(); i++)
				{
					//cout<<"i"<<i<<endl;
					float d = sensor_fusion[i][6];
					//cout<<"d"<<d<<endl;
					double vx = sensor_fusion[i][3];
					double vy = sensor_fusion[i][4];
					if(vx > 0 or vy > 0)
					{
						check_speed = sqrt(vx * vx + vy * vy);
					}
					else if(vx == 0 & vy == 0) //test to see if causing problem when car stops
					{
						check_speed = 0;
					}
					//cout<<"speed"<<check_speed<<endl;
					double check_car_s = sensor_fusion[i][5];
					//cout<<"s"<<check_car_s<<endl;
					if(d < (2 + 4 * lane + 2) & d > (2 + 4 * lane - 2)) //too close trigger
					{
						
						check_car_s += ((double)prev_size * .02 * check_speed);
						if((check_car_s > car_s) && ((check_car_s - car_s) < 35)) //was 30, 35 works well with 70 for lane changes
						{//crashing in rear on some passes at 40
						//ref_vel = 29.5; from walkthru video
							double slow_car = check_speed;
							too_close = true;
						}
					}
				
					
					if(too_close) //'or' with speed is aggressive, causing decision issues
					{
						if(lane == 0) //change left lane to middle lane
						{
							
							if(d>4 & d<8 & (check_car_s-car_s) >=0)
							{
								cost_ahead_middle = check_car_s-car_s;
								if(cost_ahead_middle < cost_benefit_middle_ahead)
								{
								
									cost_benefit_middle_ahead =  cost_ahead_middle;
								}
								else
								{
								
									cost_benefit_middle_ahead = cost_benefit_middle_ahead;
								}
								no_cars_middle_ahead = 1;
								cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
							}
							
							else if(d>4 & d<8 & (car_s-check_car_s) >=0)
							{
								cost_behind_middle = car_s-check_car_s;
								if(cost_behind_middle < cost_benefit_middle_behind)
								{
								
									cost_benefit_middle_behind = cost_behind_middle;
								}
								else
								{
								
									cost_benefit_middle_behind = cost_benefit_middle_behind;
								}
								no_cars_middle_behind = 1;
								cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
							}
							
							else if(d>0 & d<4 & (car_s-check_car_s)>=0)
							{
								cost_behind_left = car_s-check_car_s;
								if(cost_behind_left < cost_benefit_left_behind)
								{
								
									cost_benefit_left_behind = cost_behind_left;
								}
								else
								{
								
									cost_benefit_left_behind = cost_benefit_left_behind;
								}
								no_cars_left_behind = 1;
								cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
							}
								
							if(i == (sensor_fusion.size()-1)) //full sensor sweep
							{
								if(cost_benefit_middle_ahead>70 & cost_benefit_middle_behind>70 & ave_car_speed>=30)
								{
									lane = 1; //change to middle lane if clearance
									cout<<"A lane 0 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
							
								else if(cost_benefit_middle_ahead>=70 &(cost_benefit_middle_behind>=0 & cost_benefit_middle_behind<70))//had ave_car_speed<30
								{
									lane = 0; //stay in lane, not clear
								}
								
								else if((cost_benefit_middle_ahead >=0 & cost_benefit_middle_ahead<70) & cost_benefit_middle_behind>=70)
								{
									lane = 0; //stay in lane, not clear
								}
								
								else if(no_cars_middle_ahead == 0 & no_cars_middle_behind == 0 & ave_car_speed >=30)
								{
									lane = 1; //no cars ahead or behind, change to middle lane
									cout<<"B lane 0 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if(no_cars_middle_ahead == 0 & cost_benefit_middle_behind>=70 & ave_car_speed >=30)
								{
									lane = 1; //no cars ahead, clear behind, change to middle lane
									cout<<"C lane 0 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if(no_cars_middle_ahead == 0 & (cost_benefit_middle_behind>=0 & cost_benefit_middle_behind<70))
								{
									lane = 0; //stay, middle behind not clear
								}
								
								else if(cost_benefit_middle_ahead>=70 & no_cars_middle_behind == 0 & ave_car_speed >=30)
								{
									lane = 1; //front clear, no cars behind, change to middle lane
									cout<<"D lane 0 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if((cost_benefit_middle_ahead>=0 & cost_benefit_middle_ahead<70) & cost_benefit_middle_behind == 0)
								{
									lane = 0;//stay, ahead in middle lane not clear
								}
								
								else if((cost_benefit_middle_ahead>=0 & cost_benefit_middle_ahead<70)& (cost_benefit_middle_behind>=0 && cost_benefit_middle_behind<70))
								{
									lane = 0; //stay, middle ahead and behind not clear
								}
								
								else if(cost_benefit_middle_ahead < 10 || cost_benefit_middle_behind < 10)
								{
									if(d > 3 & d < 5) //crash imminent, side car lane change
									{
										ref_vel -= 0.27; //emergency braking
									}
									else
									{
										lane = 0; //stay in lane
									}
								}
								
								else if(cost_benefit_left_behind < 5 & car_speed < 5)
								{
									ref_vel +=0.20; //rear end collision imminent, speed up
								}
									
								else
								{
									lane = 0;
								}
							}
						} 
						
						else if(lane == 2) //change right lane to middle lane
						{
							
							if(d>4 & d<8 & (check_car_s-car_s) >=0)
							{
								cost_ahead_middle = check_car_s-car_s;
								if(cost_ahead_middle < cost_benefit_middle_ahead)
								{
								
									cost_benefit_middle_ahead = cost_ahead_middle;
								}
								else
								{
								
									cost_benefit_middle_ahead = cost_benefit_middle_ahead;
								}
								no_cars_middle_ahead = 1;
								cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
							}
							
							else if(d>4 & d<8 & (car_s-check_car_s) >=0)
							{
								cost_behind_middle = car_s-check_car_s;
								if(cost_behind_middle < cost_benefit_middle_behind)
								{
								
									cost_benefit_middle_behind = cost_behind_middle;
								}
								else
								{
								
									cost_benefit_middle_behind = cost_benefit_middle_behind;
								}
								no_cars_middle_behind = 1;
								cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
							}
							
							else if(d>8 & d<12 & (car_s-check_car_s)>=0)
							{
								cost_behind_right = car_s-check_car_s;
								if(cost_behind_right < cost_benefit_right_behind)
								{
								
									cost_benefit_right_behind = cost_behind_right;
								}
								else
								{
								
									cost_benefit_right_behind = cost_benefit_right_behind;
								}
								no_cars_right_behind = 1;
								cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
							}
						
							if(i == (sensor_fusion.size()-1)) //full sensor sweep
							{
								if(cost_benefit_middle_ahead>=70 & cost_benefit_middle_behind>=70 & ave_car_speed>=30)
								{
									lane = 1; //change lane to middle lane if clearance
									cout<<"A lane 2 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
							
								else if(cost_benefit_middle_ahead>=70 & (cost_benefit_middle_behind>=0 & cost_benefit_middle_behind<70))//had ave_car_speed<30
								{
									lane = 2; //too close behind for lane change, stay in lane
								}
								
								else if((cost_benefit_middle_ahead>=0 & cost_benefit_middle_ahead<70) & cost_benefit_middle_behind>=70)
								{
									lane = 2; //too close ahead for lane change, stay in lane
								}
								
								else if(no_cars_middle_ahead == 0 & no_cars_middle_behind == 0 & ave_car_speed >=30)
								{
									lane = 1; //no cars ahead or behind, change to middle lane
									cout<<"B lane 2 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if(cost_benefit_middle_ahead>=70 & no_cars_middle_behind == 0 & ave_car_speed >=30)
								{
									lane = 1; //front clear, no cars behind, change to middle lane
									cout<<"C lane 2 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if((cost_benefit_middle_ahead>=0 & cost_benefit_middle_ahead<70) & cost_benefit_middle_behind == 0)
								{
									lane = 2; //stay, not clear ahead in middle lane
								}
								
								else if(no_cars_middle_ahead == 0 & cost_benefit_middle_behind>=70 & ave_car_speed >= 30)
								{
									lane = 1; //no cars ahead, rear clear, change to middle lane
									cout<<"D lane 2 to 1"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_middle_ahead"<<cost_benefit_middle_ahead<<endl;
									cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
								}
								
								else if(no_cars_middle_ahead == 0 & (cost_benefit_middle_behind>=0 & cost_benefit_middle_behind<70))
								{
									lane = 2; //stay, not clear behind in middle lane
								}
								
								else if((cost_benefit_middle_ahead>=0 & cost_benefit_middle_ahead<70) & (cost_benefit_middle_behind>=0 && cost_benefit_middle_behind<70))
								{
									lane = 2; //stay, middle not clear ahead or behind
								}
								
								else if(cost_benefit_middle_ahead < 10 || cost_benefit_middle_behind < 10)
								{
									if(d > 7 & d < 9) //crash imminent, side car lane change
									{
										ref_vel -= 0.27; //emergency braking
									}
									else
									{
										lane = 2; //stay in lane
									}
								}
								
								else if(cost_benefit_right_behind < 5 & car_speed < 5)
								{
									ref_vel+= 0.20; //rear end collision imminent, speed up
								}
										
								else
								{
									lane = 2;
								}
							}
						}
						
						else if(lane == 1) //middle lane, change to left or right lane
						{
							if(d>0 & d<4 & (check_car_s-car_s) >=0)
							{
								cost_ahead_left = check_car_s - car_s;
								if(cost_ahead_left < cost_benefit_left_ahead)
								{
								
									cost_benefit_left_ahead = cost_ahead_left;
								}
								
								else
								{
								
									cost_benefit_left_ahead = cost_benefit_left_ahead;
								}
								no_cars_left_ahead = 1;
								cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
							}
							
							else if(d>0 & d<4 & (car_s-check_car_s) >=0)
							{
								cost_behind_left = car_s-check_car_s;
								if(cost_behind_left < cost_benefit_left_behind)
								{
								
									cost_benefit_left_behind = cost_behind_left;
								}
								else
								{
								
									cost_benefit_left_behind = cost_benefit_left_behind;
								}
								no_cars_left_behind = 1;
								cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
							}
							
							else if(d>8 & d<12 & (check_car_s-car_s) >=0)
							{
								cost_ahead_right = check_car_s-car_s;
								if(cost_ahead_right < cost_benefit_right_ahead)
								{
								
									cost_benefit_right_ahead = cost_ahead_right;
								}
								else
								{
								
									cost_benefit_right_ahead = cost_benefit_right_ahead;
								}
								no_cars_right_ahead = 1;
								cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
							}
							
							else if(d>8 & d<12 & (car_s-check_car_s) >=0)
							{
								cost_behind_right = car_s-check_car_s;
								if(cost_behind_right < cost_benefit_right_behind)
								{
								
									cost_benefit_right_behind = cost_behind_right;
								}
								else
								{
								
									cost_benefit_right_behind = cost_benefit_right_behind;
								}
								no_cars_right_behind = 1;
								cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
							}
							
							else if(d>4 & d<8 & (car_s-check_car_s) >=0)
							{
								cost_behind_middle = car_s-check_car_s;
								if(cost_behind_middle < cost_benefit_middle_behind)
								{
								
									cost_benefit_middle_behind = cost_behind_middle;
								}
								else
								{
									cost_benefit_middle_behind = cost_benefit_middle_behind;
								}
								no_cars_middle_behind = 1;
								cout<<"cost_benefit_middle_behind"<<cost_benefit_middle_behind<<endl;
							}
							
							if(i == (sensor_fusion.size()-1))//had outside lane violation at 60 and 65
							{ //full sensor sweep
								if((cost_benefit_left_ahead>=70 & cost_benefit_left_behind>=70) & (cost_benefit_right_ahead>=0 & cost_benefit_right_ahead<70)) // && car_speed > 25)
								{
									if(ave_car_speed >=30)
									{
										lane = 0; //lane change left, not clear ahead on right
										cout<<"A lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									else
									{
									
										lane = 1;
									}
									
								}
								else if(cost_benefit_left_ahead>=70 & cost_benefit_left_behind>=70 & cost_benefit_right_behind>=0 & cost_benefit_right_behind<70)
								{
									lane = 0; //lane change left, not clear behind on right
									cout<<"B lane 1 to 0"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
									cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
								}
								
								else if(cost_benefit_right_ahead>=70 & cost_benefit_right_behind>=70 & cost_benefit_left_ahead>=0 & cost_benefit_left_ahead<70) // && car_speed > 25)
								{
									lane = 2; //lane change right, left ahead not clear
									cout<<"A lane 1 to 2"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
									cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
								}
								
								else if((cost_benefit_right_ahead>=70 & cost_benefit_right_behind>=70) & (cost_benefit_left_behind>=0 & cost_benefit_left_behind<70))
								{
									lane = 2; //lane change right, left behind not clear
									cout<<"B lane 1 to 2"<<endl;
									cout<<"ave_car_speed"<<ave_car_speed<<endl;
									cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
									cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
								}
							
								else if(cost_benefit_left_ahead>=70 & cost_benefit_left_behind>=70 & cost_benefit_right_ahead>=70 & cost_benefit_right_behind>=70)
								{
									if(cost_benefit_left_ahead>=cost_benefit_right_ahead ) //failed
									{
										lane = 0;//lane change left, left more clear and usually left is fastest lane
										cout<<"C lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									else if(cost_benefit_left_ahead<cost_benefit_right_ahead & ave_car_speed >=30)
									{
										lane = 2;//lane change right, right more clear
										cout<<"C lane 1 to 2"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
										cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
									}
								}
								
								else if(no_cars_left_ahead == 0 & no_cars_right_ahead == 0)
								{//means no cars ahead except slow car
									if(cost_benefit_left_behind>=70 & cost_benefit_right_behind>=70 & ave_car_speed >=30)
									{
										lane = 0; //both lanes open, go left, usually fastest lane
										cout<<"D lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									
									else if(cost_benefit_left_behind>=70 & (cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										if(ave_car_speed >=30)
										{
											
											lane = 0; //left lane clear, go left
											cout<<"E lane 1 to 0"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
											cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
										}
										
										else
										{
										
											lane = 1;
										}
										
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70)& cost_benefit_right_behind>=70)
									{
										if(no_cars_right_behind == 1 & ave_car_speed >=30)
										{
											
											lane = 2; //right lane clear, go right
											cout<<"D lane 1 to 2"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
											cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
										}
										
										else
										{
										
											lane = 1;
										}
										
									}
									
									else if(no_cars_left_behind == 0 & ave_car_speed >=30)
									{
										
										lane = 0; //no cars front and rear, go left
										cout<<"F lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									
									else if(no_cars_right_behind == 0 & ave_car_speed >=30)
									{
										
										lane = 2; //no cars front and rear, go right
										cout<<"E lane 1 to 2"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
										cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70)&(cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										lane = 1; //stay, not clear behind left or right
									}
								}
								
								else if(no_cars_left_ahead == 0 & cost_benefit_right_ahead >= 0)
								{
									if(cost_benefit_left_behind>=70 & ave_car_speed >=30)
									{
										lane = 0; //left lane open and rear clear, go left
										cout<<"G lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									
									else if(no_cars_left_behind == 0 & ave_car_speed >=30)
									{
										
										lane = 0; //no cars behind, go left
										cout<<"H lane 1 to 0"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
										cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70) & (cost_benefit_right_ahead>=70 & cost_benefit_right_behind>=70))
									{
										if(ave_car_speed >=30)
										{
											
											lane = 2; //not clear behind on left, clear on right, go right
											cout<<"F lane 1 to 2"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
											cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
										}
										else
										{
										
											lane = 1;
										}
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70) & (cost_benefit_right_ahead>=70 & no_cars_right_behind == 0))
									{
										if(ave_car_speed >=30)
										{
											
											lane = 2; //right is clear, go right
											cout<<"G lane 1 to 2"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
											cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
										}
										
										else
										{
										
											lane = 1;
										}
										
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70)&(cost_benefit_right_ahead>=0 & cost_benefit_right_ahead<70))
									{
										lane = 1; //stay, not clear behind on left, not clear ahead on right
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70)&(cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										lane = 1; //stay, not clear behind for left or right
									}
								}
								
								else if(cost_benefit_left_ahead >= 0 & no_cars_right_ahead == 0)
								{
									if(cost_benefit_right_behind>=70 & ave_car_speed >=30)
									{
										lane = 2; //right lane open and rear clear, go right
										cout<<"H lane 1 to 2"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
										cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
									}
									
									else if(no_cars_right_behind == 0 & ave_car_speed >=30)
									{
										
										lane = 2; //no cars behind, go right
										cout<<"I lane 1 to 2"<<endl;
										cout<<"ave_car_speed"<<ave_car_speed<<endl;
										cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
										cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
									}
									
									else if((cost_benefit_left_ahead>=70 & cost_benefit_left_behind>=70) & (cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										if(ave_car_speed >=30)
										{
											
											lane = 0; //clear on left, go left, not clear right behind
											cout<<"I lane 1 to 0"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
											cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
										}
										
										else
										{
										
											lane = 1;
										}
										
									}
									
									else if((cost_benefit_left_ahead>=0 & cost_benefit_left_ahead<70) & (cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										lane = 1; //stay, not clear ahead for left or behind for right
									}
									
									else if((cost_benefit_left_behind>=0 & cost_benefit_left_behind<70)&(cost_benefit_right_behind>70))
									{
										if(ave_car_speed >=30)
										{
											lane = 2;
											cout<<"J lane 1 to 2"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_right_ahead"<<cost_benefit_right_ahead<<endl;
											cout<<"cost_benefit_right_behind"<<cost_benefit_right_behind<<endl;
										}
										else
										{
											lane = 1;
										}
									}
									
									else if((cost_benefit_left_ahead>=70 & no_cars_left_behind == 0) & (cost_benefit_right_behind>=0 & cost_benefit_right_behind<70))
									{
										if(ave_car_speed >=30)
										{
											lane = 0; //go left, left clear, right rear not clear
											cout<<"J lane 1 to 0"<<endl;
											cout<<"ave_car_speed"<<ave_car_speed<<endl;
											cout<<"cost_benefit_left_ahead"<<cost_benefit_left_ahead<<endl;
											cout<<"cost_benefit_left_behind"<<cost_benefit_left_behind<<endl;
										}
										else
										{
											lane = 1;
										}
									}
								}
									
								else if(cost_benefit_left_ahead < 10 || cost_benefit_left_behind < 10)
								{
									if(d > 3 & d < 5) //collision imminent from side car lane change
									{
										ref_vel -= 0.27; //emergency braking
									}
									else
									{
										lane = 1; //stay in lane
									}
								}
								
								else if(cost_benefit_right_ahead < 10 || cost_benefit_right_behind < 10)
								{
									if(d > 7 & d < 9) //collision imminent from side car lane change
									{
										ref_vel -= 0.27; //emergency braking
									}
									else
									{
										lane = 1; //stay in lane
									}
								}
								
								else if(cost_benefit_middle_behind < 5 & car_speed < 5)
								{
									ref_vel+= 0.20; //rear end collision imminent, speed up
								}
										
								else
								{
									lane = 1;
								}
							}
						}
					}
				}

				if(too_close)
				{ 
					if(car_speed > slow_car)
					{
						ref_vel -= 0.27;
					}
					else if(car_speed <= slow_car)
					{
						ref_vel += 0.15;
					}
					else if(car_speed <= 3.0)
					{
						ref_vel += 0.05;
					}
				}
				else if(ref_vel < 49.5)//too high at 49.75 violation
				{
					ref_vel += 0.3; //0.25 sluggish, originally 0.224
				}
				
				//set up points for path to be used with the spline, from walkthru video
				vector<double>ptsx; 
				vector<double>ptsy; 
				
				//reference state
				double ref_x = car_x;
				double ref_y = car_y;
				double ref_yaw = deg2rad(car_yaw);
			
				if (prev_size < 2)  //use car location as reference with tangent for previous pt
				{
					double prev_car_x = car_x - cos(car_yaw); //? just cos
					double prev_car_y = car_y - sin(car_yaw); //? just sin
					ptsx.push_back(prev_car_x);
					ptsx.push_back(car_x);
					ptsy.push_back(prev_car_y);
					ptsy.push_back(car_y);
				}
				
				else //more points in path, use last two pts as reference then
				{
					ref_x = previous_path_x[prev_size - 1];
					ref_y = previous_path_y[prev_size - 1];
					double ref_x_prev = previous_path_x[prev_size - 2];
					double ref_y_prev = previous_path_y[prev_size - 2];
					ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev); //angle car was heading
				
					ptsx.push_back(ref_x_prev);
					ptsx.push_back(ref_x);
					ptsy.push_back(ref_y_prev);
					ptsy.push_back(ref_y);
				}
				
				//put 3 points at 30, 60, 90 meters ahead in car lane going from Frenet to xy coordinates
				vector<double>next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
				vector<double>next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
				vector<double>next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			
				ptsx.push_back(next_wp0[0]);
				ptsx.push_back(next_wp1[0]);
				ptsx.push_back(next_wp2[0]);
			
				ptsy.push_back(next_wp0[1]);
				ptsy.push_back(next_wp1[1]);
				ptsy.push_back(next_wp2[1]);
			
				for(int i = 0; i < ptsx.size(); i++)
				{
					double shift_x = ptsx[i] - ref_x; //transformation to car reference
					double shift_y = ptsy[i] - ref_y; //transformation to car reference
				
					ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
					ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
				}
			
				tk::spline s; //create spline
			
				s.set_points(ptsx,ptsy); //set spline pts using 5 anchor pts
			
				
			
          	
				//pts for the planner
				vector<double> next_x_vals; 
				vector<double> next_y_vals;

				for(int i = 0; i < previous_path_x.size(); i++)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
			
				double target_x = 30.0;
				double target_y = s(target_x); //call spline to calculate y, given x
				double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));
			
				double x_add_on = 0;
			
				//interpolate 50 points using target distance, adding to leftover pts 
				for (int i = 1; i <= 50 - previous_path_x.size(); i++)
				{
					double N = (target_dist/(0.02 * ref_vel/2.24)); //#pts as determined by speed
					double x_point = x_add_on + (target_x)/N;
					double y_point = s(x_point); //call spline to calculate y, given x
					x_add_on = x_point;
				
					double x_ref = x_point;
					double y_ref = y_point;
				
					x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw)); //back to map reference
					y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw)); //back to map reference
				
					x_point += ref_x;
					y_point += ref_y;
				
					next_x_vals.push_back(x_point);
					next_y_vals.push_back(y_point);
				}



				// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
				json msgJson;
				msgJson["next_x"] = next_x_vals;
				msgJson["next_y"] = next_y_vals;

				auto msg = "42[\"control\","+ msgJson.dump()+"]";

				//this_thread::sleep_for(chrono::milliseconds(1000));
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			} 
        
		}
		else
		{
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
