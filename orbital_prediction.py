#!/usr/bin/env python3
"""
Open Ballistic Library - Orbital Prediction Module

This module provides functions for predicting orbital motion of objects
in space, including Keplerian orbital elements, orbital propagation,
and perturbation calculations.
"""

import math
import numpy as np
from typing import Tuple, List, Dict, Optional
from datetime import datetime, timedelta

class OrbitalPredictor:
    """
    Class for predicting orbital motion of satellites and space objects
    """
    
    def __init__(self):
        # Physical constants
        self.G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
        self.EARTH_MASS = 5.972e24  # kg
        self.EARTH_RADIUS = 6.371e6  # meters
        self.EARTH_MU = self.G * self.EARTH_MASS  # Standard gravitational parameter (m^3/s^2)
        self.J2 = 1.08263e-3  # Earth's J2 gravitational harmonic coefficient
        
    def keplerian_to_cartesian(self, semi_major_axis: float, eccentricity: float, 
                              inclination: float, raan: float, arg_perigee: float, 
                              true_anomaly: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert Keplerian orbital elements to Cartesian position and velocity vectors
        
        Args:
            semi_major_axis: Semi-major axis in meters
            eccentricity: Eccentricity (dimensionless)
            inclination: Inclination in radians
            raan: Right ascension of ascending node in radians
            arg_perigee: Argument of perigee in radians
            true_anomaly: True anomaly in radians
        
        Returns:
            Tuple of (position_vector, velocity_vector) in meters and meters/second
        """
        # Calculate position and velocity in orbital plane
        r_magnitude = semi_major_axis * (1 - eccentricity**2) / (1 + eccentricity * math.cos(true_anomaly))
        
        # Position in orbital plane
        r_x = r_magnitude * math.cos(true_anomaly)
        r_y = r_magnitude * math.sin(true_anomaly)
        
        # Velocity in orbital plane
        h = math.sqrt(self.EARTH_MU * semi_major_axis * (1 - eccentricity**2))  # Angular momentum
        v_x = self.EARTH_MU / h * (-math.sin(true_anomaly))
        v_y = self.EARTH_MU / h * (eccentricity + math.cos(true_anomaly))
        
        # Rotation matrix elements for transformation to inertial frame
        cos_raan = math.cos(raan)
        sin_raan = math.sin(raan)
        cos_inc = math.cos(inclination)
        sin_inc = math.sin(inclination)
        cos_arg = math.cos(arg_perigee)
        sin_arg = math.sin(arg_perigee)
        
        # Transform position
        pos_x = (cos_raan * cos_arg - sin_raan * sin_arg * cos_inc) * r_x - (cos_raan * sin_arg + sin_raan * cos_arg * cos_inc) * r_y
        pos_y = (sin_raan * cos_arg + cos_raan * sin_arg * cos_inc) * r_x - (sin_raan * sin_arg - cos_raan * cos_arg * cos_inc) * r_y
        pos_z = (sin_arg * sin_inc) * r_x + (cos_arg * sin_inc) * r_y
        
        # Transform velocity
        vel_x = (cos_raan * cos_arg - sin_raan * sin_arg * cos_inc) * v_x - (cos_raan * sin_arg + sin_raan * cos_arg * cos_inc) * v_y
        vel_y = (sin_raan * cos_arg + cos_raan * sin_arg * cos_inc) * v_x - (sin_raan * sin_arg - cos_raan * cos_arg * cos_inc) * v_y
        vel_z = (sin_arg * sin_inc) * v_x + (cos_arg * sin_inc) * v_y
        
        position = np.array([pos_x, pos_y, pos_z])
        velocity = np.array([vel_x, vel_y, vel_z])
        
        return position, velocity
    
    def cartesian_to_keplerian(self, position: np.ndarray, velocity: np.ndarray) -> Dict[str, float]:
        """
        Convert Cartesian position and velocity vectors to Keplerian orbital elements
        
        Args:
            position: Position vector in meters
            velocity: Velocity vector in m/s
        
        Returns:
            Dictionary containing Keplerian orbital elements
        """
        r_vec = np.array(position)
        v_vec = np.array(velocity)
        
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        
        # Specific angular momentum
        h_vec = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_vec)
        
        # Eccentricity vector
        e_vec = np.cross(v_vec, h_vec) / self.EARTH_MU - r_vec / r
        
        # Inclination
        inclination = math.acos(h_vec[2] / h)
        
        # Node vector
        n_vec = np.cross(np.array([0, 0, 1]), h_vec)
        n = np.linalg.norm(n_vec)
        
        # Right ascension of ascending node
        if n != 0:
            raan = math.acos(n_vec[0] / n)
            if n_vec[1] < 0:
                raan = 2 * math.pi - raan
        else:
            raan = 0.0
        
        # Argument of periapsis
        e = np.linalg.norm(e_vec)
        if n != 0 and e != 0:
            arg_perigee = math.acos(np.dot(n_vec, e_vec) / (n * e))
            if e_vec[2] < 0:
                arg_perigee = 2 * math.pi - arg_perigee
        else:
            arg_perigee = 0.0
        
        # True anomaly
        if e != 0:
            true_anomaly = math.acos(np.dot(e_vec, r_vec) / (e * r))
            if np.dot(r_vec, v_vec) < 0:
                true_anomaly = 2 * math.pi - true_anomaly
        else:
            true_anomaly = 0.0  # For circular orbits
        
        # Semi-major axis
        energy = v**2 / 2 - self.EARTH_MU / r
        semi_major_axis = -self.EARTH_MU / (2 * energy)
        
        return {
            'semi_major_axis': semi_major_axis,
            'eccentricity': e,
            'inclination': inclination,
            'raan': raan,  # Right ascension of ascending node
            'arg_perigee': arg_perigee,  # Argument of perigee
            'true_anomaly': true_anomaly
        }
    
    def propagate_orbit(self, position: np.ndarray, velocity: np.ndarray, 
                       time_span: float, time_step: float = 60.0) -> List[Tuple[np.ndarray, np.ndarray, float]]:
        """
        Propagate an orbit forward in time using numerical integration
        
        Args:
            position: Initial position vector in meters
            velocity: Initial velocity vector in m/s
            time_span: Total propagation time in seconds
            time_step: Time step for integration in seconds
        
        Returns:
            List of (position, velocity, time) tuples
        """
        pos = np.array(position, dtype=float)
        vel = np.array(velocity, dtype=float)
        
        results = [(pos.copy(), vel.copy(), 0.0)]
        t = 0.0
        
        while t < time_span:
            # Calculate acceleration due to gravity (with J2 perturbation)
            r = np.linalg.norm(pos)
            r_unit = pos / r
            
            # Basic gravitational acceleration
            acc = -self.EARTH_MU / r**3 * pos
            
            # J2 perturbation (simplified)
            if r > self.EARTH_RADIUS:  # Only apply if not inside Earth
                r_squared = r**2
                factor = 1.5 * self.J2 * self.EARTH_MU * (self.EARTH_RADIUS**2) / (r**5)
                z_squared = pos[2]**2
                factor_z = factor * (1 - 5 * z_squared / r_squared)
                
                acc[0] += factor_z * pos[0]
                acc[1] += factor_z * pos[1]
                acc[2] += factor * pos[2] * (3 - 5 * z_squared / r_squared)
            
            # Update velocity and position using simple integration
            vel += acc * time_step
            pos += vel * time_step
            t += time_step
            
            results.append((pos.copy(), vel.copy(), t))
        
        return results
    
    def predict_ground_track(self, position: np.ndarray, velocity: np.ndarray, 
                           time_span: float, time_step: float = 60.0) -> List[Tuple[float, float, float]]:
        """
        Predict the ground track of an orbiting object
        
        Args:
            position: Initial position vector in meters (ECI frame)
            velocity: Initial velocity vector in m/s (ECI frame)
            time_span: Total prediction time in seconds
            time_step: Time step for prediction in seconds
        
        Returns:
            List of (latitude, longitude, altitude) tuples in degrees and meters
        """
        orbit_data = self.propagate_orbit(position, velocity, time_span, time_step)
        ground_track = []
        
        # Earth rotation rate in radians per second
        earth_rate = 2 * math.pi / (24 * 3600)
        
        for pos, vel, t in orbit_data:
            x, y, z = pos
            
            # Calculate latitude
            lat = math.atan2(z, math.sqrt(x**2 + y**2))
            
            # Calculate longitude (accounting for Earth's rotation)
            lon = math.atan2(y, x) - earth_rate * t
            # Normalize longitude to [-π, π]
            while lon > math.pi:
                lon -= 2 * math.pi
            while lon < -math.pi:
                lon += 2 * math.pi
            
            # Calculate altitude
            alt = math.sqrt(x**2 + y**2 + z**2) - self.EARTH_RADIUS
            
            ground_track.append((math.degrees(lat), math.degrees(lon), alt))
        
        return ground_track
    
    def calculate_orbital_period(self, semi_major_axis: float) -> float:
        """
        Calculate the orbital period using Kepler's third law
        
        Args:
            semi_major_axis: Semi-major axis in meters
        
        Returns:
            Orbital period in seconds
        """
        return 2 * math.pi * math.sqrt(semi_major_axis**3 / self.EARTH_MU)


def predict_orbit_example():
    """
    Example function demonstrating orbital prediction capabilities
    """
    print("Open Ballistic Library - Orbital Prediction Example")
    print("=" * 50)
    
    predictor = OrbitalPredictor()
    
    # Example: ISS-like orbit
    # Semi-major axis: ~6778 km (altitude ~400 km)
    # Eccentricity: ~0.0004 (nearly circular)
    # Inclination: ~51.6 degrees
    # For this example, we'll start with a position and velocity vector
    
    # Position vector (in meters) - approximately at perigee of a 400km circular orbit
    pos = np.array([6778000.0, 0.0, 0.0])  # At equator, x-direction
    # Velocity for circular orbit at 400km altitude: ~7.67 km/s
    vel = np.array([0.0, 7670.0, 0.0])  # In y-direction (tangential)
    
    print(f"Initial position: [{pos[0]/1000:.1f}, {pos[1]/1000:.1f}, {pos[2]/1000:.1f}] km")
    print(f"Initial velocity: [{vel[0]:.1f}, {vel[1]:.1f}, {vel[2]:.1f}] m/s")
    
    # Convert to Keplerian elements to verify
    elements = predictor.cartesian_to_keplerian(pos, vel)
    print(f"\nDerived Keplerian elements:")
    print(f"  Semi-major axis: {elements['semi_major_axis']/1000:.1f} km")
    print(f"  Eccentricity: {elements['eccentricity']:.4f}")
    print(f"  Inclination: {math.degrees(elements['inclination']):.2f}°")
    print(f"  Orbital period: {predictor.calculate_orbital_period(elements['semi_major_axis'])/60:.2f} min")
    
    # Propagate orbit for 1 period
    period = predictor.calculate_orbital_period(elements['semi_major_axis'])
    orbit_path = predictor.propagate_orbit(pos, vel, period, time_step=120.0)
    
    print(f"\nOrbit propagated for {len(orbit_path)} time steps")
    print(f"Predicted orbital period: {period:.2f} seconds ({period/60:.2f} minutes)")
    
    # Calculate ground track for 1/4 of the orbit
    ground_track = predictor.predict_ground_track(pos, vel, period/4, time_step=60.0)
    print(f"\nGround track points for 1/4 orbit: {len(ground_track)}")
    if ground_track:
        start_point = ground_track[0]
        end_point = ground_track[-1]
        print(f"Start: Lat {start_point[0]:.2f}°, Lon {start_point[1]:.2f}°")
        print(f"End:   Lat {end_point[0]:.2f}°, Lon {end_point[1]:.2f}°")


if __name__ == "__main__":
    predict_orbit_example()