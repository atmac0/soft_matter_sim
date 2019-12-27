#include "collisions.h"

/* PRECONDITION  : None */
/* POSTCONDITION : A collision object will have been initialized. The collisions array will have been sized accordingly depending on particle type. */
/* */
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

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* Returns the number of collisions in the list. */
uint32_t  Collisions::counts()
{
  return collision_counts;
}

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* */
collision_t Collisions::get_collision_at(uint32_t index)
{
  return collisions[index];
}

/* PRECONDITION  : None */
/* POSTCONDITION : A collisions will have been added to the collisions array. The collisions counter will be incremented. */
/* */
void Collisions::add_collision(collision_t collision)
{
  collisions[collision_counts] = collision;
  collision_counts++;
}

/* PRECONDITION  : None */
/* POSTCONDITION : The collisions array will have been freed. */
/* */
Collisions::~Collisions()
{
  free(collisions);
}
