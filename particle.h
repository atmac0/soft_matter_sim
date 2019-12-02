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
 
  uint16_t particle_num;            //unique number to use as a particle identifier
  Particle_type particle_type;      //type of particle, used to identify shape.

  double angular_velocity;          //current angular velocity, described in rads/micro_second
  double orientation;               //orientation measured as an angle in radians. Describes the current angle of the object,
                                    //relative to the non-rotated orientation.

  double x_velocity;                //translational velocity in the x direction, measured as left to right as the increasing direction, in micro-meters per micro-second
  double y_velocity;                //translational velocity in the y direction, measured as up to down as the increasing direction in micro-meters per micro-second

  uint32_t moment_of_inertia;       //moment of inertia. dependent on mass and size of particle

  coord_t center_mass_coord;          //coordinate in the field of the center of mass
  double x_position_granular;      //x position within a cell, setting the origin at the left side of the cell
  double y_position_granular;      //y position within a cell, setting the origin at the top of the cell
  
  Field* field;
  
  std::vector<coord_t> edge_locations;           //list of all the coordinates of the location of edges. To be used for the deletion of edged
  
public:
  Particle(Field* field_, uint16_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord);
  int32_t hypotenuse(int32_t x, int32_t y);
  Collisions* draw_square();
  Collisions* draw_edges();
  void delete_edges();
  uint16_t get_particle_mass();
  void calc_x_y_from_angle(int32_t &x, int32_t &y, double angle);
  coord_t translate_to_field(coord_t edge_point);
  void add_edge_location(coord_t);
  void propagate();
  Collisions* translate_x_by_1();
  Collisions* translate_y_by_1();
  void resolve_collisions(Collisions* collisions);
  uint32_t translate_y_by_granular(double granularity);
  uint32_t translate_x_by_granular(double granularity);
  void rotate_particle(double time_span);
};
