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
  double orientation;               //orientation measured as an angle in radians. Describes the current angle of the object,
                                    //relative to the non-rotated orientation.

  uint32_t trans_velocity;          //translational velocity measured in cells per cycle
  double direction;                 //translational direction measured as an angle in radians

  double angular_velocity;          //current angular velocity, described in rads/second
  uint32_t moment_of_inertia;       //moment of inertia. dependent on mass and size of particle

  coord_t center_mass_coord;          //coordinate in the field of the center of mass

  Field* field;

  coord_t edge_locations[];           //list of all the coordinates of the location of edges. To be used for the deletion of edged

public:
  Particle(Field* field_, uint16_t par_num, Particle_type par_type, double orient, uint32_t trans_vel, double angular_vel, double dir, coord_t cm_coord);
  int32_t hypotenuse(int32_t x, int32_t y);
  Collisions draw_square();
  Collisions draw_edges();
  void delete_edges();
  uint16_t get_particle_mass();
  void calc_x_y_from_angle(int32_t &x, int32_t &y, double angle);
  coord_t translate_to_field(coord_t edge_point);
};
