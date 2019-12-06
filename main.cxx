#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "data_types.h"
#include <string>
#include <memory>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "field.h"
#include "particle.h"

int32_t main()
{
 
  Field my_field;
  
  field_t par_num0  = 0;
  Particle_type par0 = SQUARE;
  double ang_vel0    = 0;
  double orient0     = 0;
  double x_vel0      = 10;
  double y_vel0      = 0;
  coord_t cm_coord0  = {300,500};

  field_t par_num1  = 1;
  Particle_type par1 = SQUARE;
  double ang_vel1    = 0;
  double orient1     = 0;
  double x_vel1      = -10;
  double y_vel1      = 0;
  coord_t cm_coord1 = {500,500};  

  field_t par_num2  = 2;
  Particle_type par2 = SQUARE;
  double ang_vel2    = -.25;
  double orient2     = M_PI/3;
  double x_vel2      = 27;
  double y_vel2      = 0;
  coord_t cm_coord2 = {500,700};  
  
  Particle my_particle0(&my_field, par_num0, par0, ang_vel0, orient0, x_vel0, y_vel0, cm_coord0);
  Particle my_particle1(&my_field, par_num1, par1, ang_vel1, orient1, x_vel1, y_vel1, cm_coord1); 
  Particle my_particle2(&my_field, par_num2, par2, ang_vel2, orient2, x_vel2, y_vel2, cm_coord2);

  Particle* particles[NUM_PARTICLES];
  particles[0] = &my_particle0;
  particles[1] = &my_particle1;
  //particles[2] = &my_particle2;

  for(int i = 0; i < NUM_PARTICLES; i++)
  {
    particles[i]->set_particles_array(particles);
    particles[i]->draw_edges();
  }
  
  for(uint32_t frame_num = 1; frame_num < 400; frame_num++)
  {
    for(int i = 0; i < NUM_PARTICLES; i++)
    {
      particles[i]->propagate();
    }
    my_field.field_to_png();
  }
  
  return 0;
}

