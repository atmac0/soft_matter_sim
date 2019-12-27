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

field_t find_lowest_time_particle(Particle* particles[NUM_PARTICLES]);
Particle make_random_square(field_t par_num, coord_t cm_coord, Field* my_field);
void initialize_particles(Particle* particles[NUM_PARTICLES], Field* my_field);

int main()
{
  std::cout <<"here";
  std::cout << std::rand();
  Field my_field;

  uint32_t frame_counter = 0;
<<<<<<< HEAD
  double time_increment = .1;
  double time = 0;
  double end_time = 1000000;
=======
  uint32_t time_increment = 5;
  uint32_t time = 0;
  uint32_t end_time = 1000000;
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.

  Particle* particles[NUM_PARTICLES];
  initialize_particles(particles, &my_field);
  
  field_t lowest_time_particle;
  
  while(time < end_time)
  {
    std::cout << "TIME: " << time << "\n";
    for(int i = 0; i < NUM_PARTICLES; i++)
    {
      lowest_time_particle = find_lowest_time_particle(particles);
      particles[lowest_time_particle]->propagate();
    }

    if(time <= particles[lowest_time_particle]->get_relative_time())
    {
      my_field.field_to_png();
      time += time_increment;
    }
  }

  return 0;
}


field_t find_lowest_time_particle(Particle* particles[NUM_PARTICLES])
{
  double time = particles[0]->get_relative_time();
  field_t par_num = 0;
  
  for(uint32_t i = 1; i < NUM_PARTICLES; i++)
  {
    if(particles[i]->get_relative_time() < time)
    {
      time = particles[i]->get_relative_time();
      par_num = i;
    }
  }
  
  return par_num;
}

Particle make_random_square(field_t par_num, coord_t cm_coord, Field* my_field)
{

  Particle_type par_type = SQUARE;
  double ang_vel    = (std::rand()/double(RAND_MAX))/15;
  double orient     = std::rand()/double(RAND_MAX);
  double x_vel      = (std::rand()/double(RAND_MAX))*5;
  double y_vel      = (std::rand()/double(RAND_MAX))*5;
  
  Particle particle(my_field, par_num, par_type, ang_vel, orient, x_vel, y_vel, cm_coord);
  return particle;
}

void initialize_particles(Particle* particles[NUM_PARTICLES], Field* my_field)
{

  uint32_t x_pos = 10;
  uint32_t y_pos = 10;
  coord_t cm_coord;
  Particle particle;
  
  for(uint32_t par_num = 0; par_num < NUM_PARTICLES; par_num++)
  {

    cm_coord.x = x_pos;
    cm_coord.y = y_pos;

    Particle_type par_type = SQUARE;
    double ang_vel    = (std::rand()/double(RAND_MAX))/15;
    double orient     = std::rand()/double(RAND_MAX);
    double x_vel      = (std::rand()/double(RAND_MAX))*5;
    double y_vel      = (std::rand()/double(RAND_MAX))*5;
    
    particle = new Particle(my_field, par_num, par_type, ang_vel, orient, x_vel, y_vel, cm_coord);
    
    particles[par_num] = &particle;

    x_pos += 40;
    if(x_pos > (FIELD_WIDTH - 20))
    {
      x_pos = 10;
      y_pos += 40;
    }
  }

  for(int i = 0; i < NUM_PARTICLES; i++)
  {
<<<<<<< HEAD
    (*it)->set_particles_vector(particles);
    (*it)->draw();
=======
    particles[i]->set_particles_array(particles);
    particles[i]->draw_edges();
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.
  }

  my_field->field_to_png();
  
}
