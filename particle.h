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

  double angular_velocity;          /* current angular velocity, described in rads/micro_second, in the clockwise direction. */
  double orientation;               /* orientation measured as an angle in radians. Describes the current angle of the object,
				      in the clockwise direction, relative to the non-rotated orientation. */

  double x_velocity;                /* translational velocity in the x direction, measured as left to right as the increasing direction,
				      in micro-meters per micro-second */
  double y_velocity;                /* translational velocity in the y direction, measured as up to down as the increasing direction in 
				       micro-meters per micro-second */

  uint32_t mass;                    /* mass in grams */
  double moment_of_inertia;         /* moment of inertia. dependent on mass and size of particle */
  
  coord_t center_mass_coord;        /* coordinate in the field of the center of mass */
  double x_position_granular;       /* x position within a cell, setting the origin at the left side of the cell */
  double y_position_granular;       /* y position within a cell, setting the origin at the top of the cell */

  double relative_time;             /* total time the particle has spent propagating */
  
  Field* field;
  
<<<<<<< HEAD
  std::vector<coord_t> edge_locations; /* vector of all the coordinates of the location of edges. To be used for the deletion of edged */
  std::vector<Particle*> particles;    /* vector shared between all particles, of all the particles in the field. */
=======
  std::vector<coord_t> edge_locations;           //list of all the coordinates of the location of edges. To be used for the deletion of edged
  Particle** particles;
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.
  
public:
  Particle(Field* field_, field_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord);
  void set_particles_vector(std::vector<Particle*>* particles);
  
  Collisions* draw();
  Collisions* draw_square();
  void delete_edges();
  
  coord_t translate_to_field(coord_t edge_point);
  double find_angle_relative_to_field(double x, double y);
  void add_edge_location(coord_t);
<<<<<<< HEAD

  void translate_x_by_1();
  void translate_y_by_1();
  void translate_y_by_granular(double granularity);
  void translate_x_by_granular(double granularity);
  void rotate_particle(double time_span);  

  void propagate();
  void resolve_collisions(Collisions* collisions);

  momentum_t find_linear_momentum_at(coord_t point, field_t particle_num);
  
=======
  void propagate();
  Collisions* translate_x_by_1();
  Collisions* translate_y_by_1();
  Collisions* rotate_particle(double time_span);
  void resolve_collisions(Collisions* collisions);
  uint32_t translate_y_by_granular(double granularity);
  uint32_t translate_x_by_granular(double granularity);
  momentum_t find_linear_momentum_at(coord_t point);
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.
  uint32_t get_mass();
  coord_t get_center_mass_coord();
  double get_x_velocity();
  double get_y_velocity();
  double get_angular_velocity();
<<<<<<< HEAD
=======
  void set_particles_array(Particle** par_arr);
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.
  double get_moment_of_inertia();
  double get_relative_time();
<<<<<<< HEAD
  field_t get_particle_num();
  double find_r_cell_speed();
  
  void increment_y_velocity(double dy);
  void increment_x_velocity(double dx);
  void increment_angular_velocity(double dw);

  bool is_change_in_linear_in_direction_of_cm(coord_t coll_location, double d_linear_x, double d_linear_y);
  bool is_change_in_angular_in_direction_of_cm(coord_t coll_location, double d_linear_x, double d_linear_y);
  int32_t determine_direction_of_angular_change(int32_t x_diff, int32_t y_diff, double dp_lx, double dp_ly);
  void find_change_in_momentum(uint16_t par1, uint16_t par2, momentum_t linear_momentum, coord_t coll_location, momentum_t dp_t[], double dp_l[]);

=======
  int32_t determine_direction_of_angular_change(int32_t x_diff, int32_t y_diff, double dp_lx, double dp_ly);
  double find_r_cell_speed();
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.
};
