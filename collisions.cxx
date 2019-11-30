#include "collisions.h"

Collisions::Collisions(Particle_type particle_type)
{

  uint32_t max_collisions;
  
  if(particle_type == SQUARE)
  {
    max_collisions = NUM_SIDES_SQUARE*(SQUARE_SIDE_LENGTH/CELL_SIZE);
  }

  collisions = (collision_t*)malloc(sizeof(collision_t)*max_collisions);
  collision_count = 0;
}

uint32_t  Collisions::counts()
{
  return collision_count;
}

collision_t Collisions::get_collision_at(uint32_t index)
{
  return collisions[index];
}

uint32_t Collisions::get_collision_count()
{
  return collision_count;
}

void Collisions::add_collision(collision_t collision)
{
  collisions[collision_count] = collision;
  collision_count++;
}

Collisions::~Collisions()
{
  free(collisions);
}
