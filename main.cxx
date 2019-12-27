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

field_t find_lowest_time_particle(std::vector<Particle*> particles);
Particle make_random_square(field_t par_num, coord_t cm_coord, Field* my_field);
void initialize_particles(std::vector<Particle*>* particles, Field* my_field);

int main()
{
  Field my_field;

  uint32_t frame_counter = 0;
  double time_increment = .1;
  double time = 0;
  double end_time = 1000000;

  std::vector<Particle*> particles;
  initialize_particles(&particles, &my_field);
  
  field_t lowest_time_particle;
  
  while(time < end_time)
  {
    std::cout << "TIME: " << time << "\n";
    for(int i = 0; i < NUM_PARTICLES; i++)
    {
      lowest_time_particle = find_lowest_time_particle(particles);
      particles.at(lowest_time_particle)->propagate();
    }

    if(time <= particles.at(lowest_time_particle)->get_relative_time())
    {
      my_field.field_to_png();
      time += time_increment;
    }
  }

  return 0;
}


field_t find_lowest_time_particle(std::vector<Particle*> particles)
{
  double time = particles.front()->get_relative_time();
  field_t par_num = 0;
  
  for(std::vector<Particle*>::iterator it = particles.begin(); it != particles.end(); it++)
  {
    if((*it)->get_relative_time() < time)
    {
      time = (*it)->get_relative_time();
      par_num = (*it)->get_particle_num();
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

void initialize_particles(std::vector<Particle*>* particles, Field* my_field)
{

  uint32_t x_pos = 10;
  uint32_t y_pos = 10;
  coord_t cm_coord;
  
  for(uint32_t par_num = 0; par_num < NUM_PARTICLES; par_num++)
  {

    cm_coord.x = x_pos;
    cm_coord.y = y_pos;

    Particle_type par_type = SQUARE;
    double ang_vel    = (std::rand()/double(RAND_MAX))/30;
    double orient     = std::rand()/double(RAND_MAX);
    double x_vel      = (std::rand()/double(RAND_MAX))*10;
    double y_vel      = (std::rand()/double(RAND_MAX))*10;
    
    particles->push_back(new Particle(my_field, par_num, par_type, ang_vel, orient, x_vel, y_vel, cm_coord));

    x_pos += 150;
    if(x_pos > (FIELD_WIDTH - 150))
    {
      x_pos = 150;
      y_pos += 150;
    }
  }

  for(std::vector<Particle*>::iterator it = particles->begin(); it != particles->end(); it++)
  {
    (*it)->set_particles_vector(particles);
    (*it)->draw();
  }

  my_field->field_to_png();
  
}
