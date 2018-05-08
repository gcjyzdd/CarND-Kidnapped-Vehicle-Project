/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// set number of particles
	num_particles = 600;
	is_initialized = true;
	default_random_engine gen;
	// create normal distribution for x, y, and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_t(theta, std[2]);

	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_t(gen);
		p.weight = 1.;
		particles.push_back(p);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	// set standard deviation of noises
	double dx = 0.2;
	double dy = 0.2;
	double dt = 0.001;
	normal_distribution<double> dist_x(0, dx);
	normal_distribution<double> dist_y(0, dy);
	normal_distribution<double> dist_t(0, dt);

	// predict theta
	for (int i = 0; i < num_particles; i++)
	{
		particles[i].theta += yaw_rate * delta_t + dist_t(gen);
	}
	// predict x and y
	if (std::abs(yaw_rate) > 1e-6)
	{
		for (int i = 0; i < num_particles; i++)
		{
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y += velocity / yaw_rate * (-cos(particles[i].theta + yaw_rate * delta_t) + cos(particles[i].theta)) + +dist_y(gen);
		}
	}
	else // yaw rate is zero
	{
		for (int i = 0; i < num_particles; i++)
		{
			particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
			particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
		}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> &predicted, std::vector<LandmarkObs> &observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	// predicted is predicted land mark positions of a specific particles
	// observations is the actual observation of landmarks
	
	double b;	// best distance
	double tmp; // teporary variable to store current distance

	for (size_t j = 0; j < predicted.size(); j++)	
	{
		b = 1e12;
		for (size_t i = 0; i < observations.size(); i++)	
		{
			// use squared distance. There is no need to apply root square here.
			tmp = pow(predicted[j].x - observations[i].x, 2) + pow(predicted[j].y - observations[i].y, 2);
			if (tmp < b)
			{// a better match is found
				b = tmp;
				// correlate this prediction to landmark i.
				predicted[j].id = i;
			}
		}
	}
}

/**
 * calculate a multi-variable distribution probability
 * @param x 	predcited x coordinate of a landmark
 * @param y 	predcited y coordinate of a landmark
 * @param mux	measured x coordinate of a landmark
 * @param muy	measured y coordinate of a landmark
 * @param sigx  standard deviation of sensor measurement regarding to x
 * @param sigy  standard deviation of sensor measurement regarding to y
 * */
double Pxy(double x, double y, double mux, double muy, double sigx, double sigy)
{
	return 1./(2. * M_1_PI * sigx * sigy) * exp( -0.5 * ( pow((x-mux)/sigx, 2) + pow((y-muy)/sigy, 2) ) );
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
								   const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	std::vector<LandmarkObs> predicted;
	// get ground truth
	std::vector<LandmarkObs> gt;
	LandmarkObs lm;
	for(size_t i=0;i<map_landmarks.landmark_list.size();i++)
	{
		lm.id = map_landmarks.landmark_list[i].id_i;
		lm.x = map_landmarks.landmark_list[i].x_f;
		lm.y = map_landmarks.landmark_list[i].y_f;
		gt.push_back(lm);
	}
	
	double t, w_tmp;

	for(int i=0;i<num_particles;i++)
	{
		predicted.clear();
		for(size_t j=0; j<observations.size();j++)
		{
			t = particles[i].theta;
			lm.x = cos(t)*observations[j].x - sin(t)*observations[j].y + particles[i].x;
			lm.y = sin(t)*observations[j].x + cos(t)*observations[j].y + particles[i].y;
			predicted.push_back(lm);
		}
		dataAssociation(predicted, gt);
		w_tmp = 1.;
		for(size_t k=0;k<predicted.size();k++)
		{
			int ind = predicted[k].id;
			// return the product of all landmarks' probabilities
			w_tmp *= Pxy(predicted[k].x, predicted[k].y, gt[ind].x, gt[ind].y, std_landmark[0], std_landmark[1]);
		}
		particles[i].weight = w_tmp;
	}
	// as long as we resample based on proportional weights, there is no need to normalize.
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
    std::mt19937 gen(rd());
	std::vector<double> w(num_particles);
	for(int i=0;i<num_particles;i++)
	{
		w[i] = particles[i].weight;
	}
	// generate a discrete distribution based on weights
    std::discrete_distribution<int> d(w.begin(), w.end());

	std::vector<Particle> tmp = particles;
	int ind;
	for(int i=0; i<num_particles;i++)
	{
		ind = d(gen);
		particles[i] = tmp[ind];
	}	

}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
										 const std::vector<double> &sense_x, const std::vector<double> &sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
