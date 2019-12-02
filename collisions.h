#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "data_types.h"

struct collision_t
{
  coord_t location;
  int16_t particle1;
  int16_t particle2;
};

class Collisions
{
 private:
  uint32_t collision_count;
  collision_t* collisions;
 public:
  Collisions(Particle_type particle_type);
  uint32_t counts();
  collision_t get_collision_at(uint32_t index);
  void add_collision(collision_t collision);
  ~Collisions();
};
