function [value,isterminal,direction] = PoincareSection(y,prey_val,stop)
  value = y(1) - prey_val; % The value that we want to be zero
  isterminal = stop;  % Halt integration 1:  continue 0
  direction = 1;   
end