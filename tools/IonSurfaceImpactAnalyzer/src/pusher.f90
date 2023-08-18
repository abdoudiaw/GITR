module pusher
    use share 
    use utilities
    implicit none

    contains

    subroutine borisPusher()
        use share
        implicit none
        
        do j = 1, npart
            ! Calculate electric field for particle j
            !!!!! write a function to get the electric field
            ex(j) = 0.0
            ey(j) = 0.0
            ez(j) = sheath_field(abs(z(j)), alpha, larmor_over_debye_radius) * zbar_imp * tempe_back / tempi_back
            ! Calculate magnetic field for particle j
            bx(j) = Bmagnitude * cos( alpha * pi / 180.0) 
            by(j) = 0.0
            bz(j) = Bmagnitude * sin( alpha * pi / 180.0)
            
            q_over_m =  1.0  ! this is one because we are using normalized units
            ! Half-step electric field update
            v_minus(1) = vx(j) + 0.5 * q_over_m * ex(j) * dt
            v_minus(2) = vy(j) + 0.5 * q_over_m * ey(j) * dt
            v_minus(3) = vz(j) + 0.5 * q_over_m * ez(j) * dt
        
            ! Magnetic rotation
            t(1) = q_over_m * bx(j) * 0.5 * dt
            t(2) = q_over_m * by(j) * 0.5 * dt
            t(3) = q_over_m * bz(j) * 0.5 * dt
              
            s(1) = 2*t(1) / (1 + dot_product(t, t))
            s(2) = 2*t(2) / (1 + dot_product(t, t))
            s(3) = 2*t(3) / (1 + dot_product(t, t))
              
            ! Compute cross product v_minus x t
            cross_v_minus_t(1) = v_minus(2) * t(3) - v_minus(3) * t(2)
            cross_v_minus_t(2) = v_minus(3) * t(1) - v_minus(1) * t(3)
            cross_v_minus_t(3) = v_minus(1) * t(2) - v_minus(2) * t(1)
   
            ! Compute v_prime
            v_prime(1) = v_minus(1) + cross_v_minus_t(1)
            v_prime(2) = v_minus(2) + cross_v_minus_t(2)
            v_prime(3) = v_minus(3) + cross_v_minus_t(3)
   
            ! Cross product calculations go here for v_prime
              
            v_plus(1) = v_minus(1) + cross_product_component(v_prime, s, 1)
            v_plus(2) = v_minus(2) + cross_product_component(v_prime, s, 2)
            v_plus(3) = v_minus(3) + cross_product_component(v_prime, s, 3)
              
            ! Another half-step electric field update for particle j
            vx(j) = v_plus(1) + 0.5 * q_over_m * ex(j) * dt
            vy(j) = v_plus(2) + 0.5 * q_over_m * ey(j) * dt
            vz(j) = v_plus(3) + 0.5 * q_over_m * ez(j) * dt
              
            ! Update positions for particle j
            x(j) = x(j) + vx(j) * dt
            y(j) = y(j) + vy(j) * dt
            z(j) = z(j) + vz(j) * dt
        end do
    end subroutine borisPusher
    
end module pusher