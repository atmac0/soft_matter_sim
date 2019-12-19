#include "collisions.h"

Collisions::Collisions(Particle_type particle_type)
{
  uint32_t max_collisions;
  
  if(particle_type == SQUARE)
  {
    max_collisions = NUM_SIDES_SQUARE*(SQUARE_SIDE_LENGTH/CELL_SIZE)*LINE_THICKNESS;
  }

  collisions = (collision_t*)malloc(sizeof(collision_t)*max_collisions);
  collision_counts = 0;
}

uint32_t  Collisions::counts()
{
  return collision_counts;
}

collision_t Collisions::get_collision_at(uint32_t index)
{
  return collisions[index];
}

void Collisions::add_collision(collision_t collision)
{
  collisions[collision_counts] = collision;
  collision_counts++;
}

Collisions::~Collisions()
{
  free(collisions);
}
