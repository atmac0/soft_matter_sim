#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "data_types.h"
#include <memory>
#include <cmath>
#include "collisions.h"
#include "field.h"

class Particle
{
private:
 
  field_t particle_num;            //unique number to use as a particle identifier
  Particle_type particle_type;      //type of particle, used to identify shape.

  double angular_velocity;          //current angular velocity, described in rads/micro_second, in the clockwise direction.
  double orientation;               /*orientation measured as an angle in radians. Describes the current angle of the object,
				      in the clockwise direction, relative to the non-rotated orientation. */

  double x_velocity;                //translational velocity in the x direction, measured as left to right as the increasing direction, in micro-meters per micro-second
  double y_velocity;                //translational velocity in the y direction, measured as up to down as the increasing direction in micro-meters per micro-second

  uint32_t mass;
  double moment_of_inertia;       //moment of inertia. dependent on mass and size of particle
  
  coord_t center_mass_coord;          //coordinate in the field of the center of mass
  double x_position_granular;      //x position within a cell, setting the origin at the left side of the cell
  double y_position_granular;      //y position within a cell, setting the origin at the top of the cell

  double relative_time;
  
  Field* field;
  
  std::vector<coord_t> edge_locations;           //list of all the coordinates of the location of edges. To be used for the deletion of edged
  std::vector<Particle*> particles;
  
public:
  Particle(Field* field_, field_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord);
  int32_t hypotenuse(int32_t x, int32_t y);
  Collisions* draw_square();
  Collisions* draw_edges();
  void delete_edges();
  uint16_t get_particle_mass();
  void calc_x_y_from_angle(int32_t &x, int32_t &y, double angle);
  coord_t translate_to_field(coord_t edge_point);
  void add_edge_location(coord_t);
  void propagate();
  void translate_x_by_1();
  void translate_y_by_1();
  void rotate_particle(double time_span);
  void resolve_collisions(Collisions* collisions);
  void translate_y_by_granular(double granularity);
  void translate_x_by_granular(double granularity);
  momentum_t find_linear_momentum_at(coord_t point);
  uint32_t get_mass();
  coord_t get_center_mass_coord();
  double get_x_velocity();
  double get_y_velocity();
  void increment_y_velocity(double dy);
  void increment_x_velocity(double dx);
  void find_change_in_momentum(uint16_t par1, uint16_t par2, momentum_t linear_momentum, coord_t cm_coord, coord_t coll_location, momentum_t dp_t[], double dp_l[]);
  momentum_t find_linear_momentum_at(coord_t point, field_t particle_num);
  double get_angular_velocity();
  void set_particles_vector(std::vector<Particle*>* particles);
  double get_moment_of_inertia();
  void increment_angular_velocity(double dw);
  double get_relative_time();
  int32_t determine_direction_of_angular_change(int32_t x_diff, int32_t y_diff, double dp_lx, double dp_ly);
  double find_r_cell_speed();
  field_t get_particle_num();
  bool is_change_in_linear_in_direction_of_cm(coord_t cm_coord, coord_t coll_location, double d_linear_x, double d_linear_y);
};
