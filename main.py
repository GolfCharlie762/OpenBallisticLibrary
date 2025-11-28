#!/usr/bin/env python3
"""
Open Ballistic Library - Main Module

This is the main entry point for the open-source ballistics library.
The library provides tools for ballistic calculations, trajectory predictions,
and orbital mechanics for projectiles and space objects.
"""

import math
import numpy as np
from typing import Tuple, List, Optional

class BallisticCalculator:
    """
    Main class for ballistic calculations including trajectory prediction,
    atmospheric effects, and orbital mechanics.
    """
    
    def __init__(self):
        # Physical constants
        self.G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
        self.EARTH_MASS = 5.972e24  # kg
        self.EARTH_RADIUS = 6.371e6  # meters
        self.EARTH_MU = self.G * self.EARTH_MASS  # Standard gravitational parameter
        
    def calculate_trajectory(self, initial_velocity: float, launch_angle: float, 
                           mass: float, drag_coefficient: float = 0.47, 
                           cross_section: float = 0.01, time_step: float = 0.01) -> List[Tuple[float, float, float]]:
        """
        Calculate projectile trajectory with atmospheric drag
        
        Args:
            initial_velocity: Initial velocity in m/s
            launch_angle: Launch angle in degrees
            mass: Mass of projectile in kg
            drag_coefficient: Drag coefficient (default for sphere)
            cross_section: Cross-sectional area in m^2
            time_step: Time step for simulation in seconds
        
        Returns:
            List of (time, x, y) positions
        """
        # Convert angle to radians
        angle_rad = math.radians(launch_angle)
        
        # Initial conditions
        vx = initial_velocity * math.cos(angle_rad)
        vy = initial_velocity * math.sin(angle_rad)
        x, y = 0.0, 0.0
        t = 0.0
        
        positions = [(t, x, y)]
        
        # Atmospheric density at sea level (kg/m^3)
        rho = 1.225
        
        while y >= 0:  # Continue until projectile hits ground
            # Calculate drag force
            v_total = math.sqrt(vx**2 + vy**2)
            drag_force = 0.5 * rho * v_total**2 * drag_coefficient * cross_section
            
            # Calculate acceleration components
            ax = -drag_force * vx / v_total / mass if v_total != 0 else 0
            ay = -drag_force * vy / v_total / mass - 9.81  # gravity
            
            # Update velocities
            vx += ax * time_step
            vy += ay * time_step
            
            # Update positions
            x += vx * time_step
            y += vy * time_step
            t += time_step
            
            positions.append((t, x, y))
        
        return positions
    
    def orbital_velocity(self, altitude: float) -> float:
        """
        Calculate circular orbital velocity at given altitude
        
        Args:
            altitude: Altitude above Earth's surface in meters
        
        Returns:
            Orbital velocity in m/s
        """
        r = self.EARTH_RADIUS + altitude
        return math.sqrt(self.EARTH_MU / r)
    
    def escape_velocity(self, altitude: float) -> float:
        """
        Calculate escape velocity at given altitude
        
        Args:
            altitude: Altitude above Earth's surface in meters
        
        Returns:
            Escape velocity in m/s
        """
        r = self.EARTH_RADIUS + altitude
        return math.sqrt(2 * self.EARTH_MU / r)


def main():
    """
    Main function demonstrating the ballistics library capabilities
    """
    print("Open Ballistic Library - Main Module")
    print("=" * 40)
    
    # Create calculator instance
    calc = BallisticCalculator()
    
    # Example: Calculate trajectory
    print("\nExample: Projectile trajectory")
    trajectory = calc.calculate_trajectory(
        initial_velocity=100.0,  # m/s
        launch_angle=45.0,      # degrees
        mass=1.0,               # kg
        drag_coefficient=0.47,  # sphere
        cross_section=0.01      # m^2
    )
    
    # Print first few and last few points
    print(f"Initial conditions: v=100 m/s, angle=45°")
    print(f"Trajectory points: {len(trajectory)}")
    print(f"Flight time: {trajectory[-1][0]:.2f} seconds")
    print(f"Range: {trajectory[-1][1]:.2f} meters")
    
    # Example: Orbital calculations
    print("\nExample: Orbital mechanics")
    altitude = 400000  # 400 km altitude (ISS-like)
    orb_vel = calc.orbital_velocity(altitude)
    esc_vel = calc.escape_velocity(altitude)
    
    print(f"At {altitude/1000:.0f} km altitude:")
    print(f"  Circular orbital velocity: {orb_vel:.2f} m/s")
    print(f"  Escape velocity: {esc_vel:.2f} m/s")


if __name__ == "__main__":
    main()