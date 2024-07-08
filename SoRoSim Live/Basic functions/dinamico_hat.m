function se3 = dinamico_hat(screw) % optimized on 30.05.2022
se3 = [0 -screw(3) screw(2) screw(4);screw(3) 0 -screw(1) screw(5);-screw(2) screw(1) 0 screw(6);0 0 0 0];
