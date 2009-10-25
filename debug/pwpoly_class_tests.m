function[container] = pwpoly_class_tests()
% [container] = pwpoly_class_tests()
%
%     Runs various tests for the matlab-implemented PiecewisePolynomial class.

global packages
import debug.*

container = TestContainer();

opt.K = 10;
