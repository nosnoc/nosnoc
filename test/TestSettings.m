classdef TestSettings < matlab.unittest.TestCase
    %TESTSETTINGS Summary of this class goes here
    %   Detailed explanation goes here    
    properties (TestParameter)
        type = {'single','double','uint16'};
        level = struct('small',2,'medium',4,'large',6);
        side = struct('small',9,'medium',81,'large',729);
    end
    
    methods (Test)
        function testRemainPixels(testCase,level)
            expPixelCount = 8^level;
            actPixels = find(sierpinski(level));
            testCase.verifyNumElements(actPixels,expPixelCount)
        end
        
        function testClass(testCase,type,level)
            testCase.verifyClass( ...
                sierpinski(level,type),type)
        end
        
        function testDefaultL1Output(testCase)
            exp = single([1 1 1; 1 0 1; 1 1 1]);
            testCase.verifyEqual(sierpinski(1),exp)
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        function testNumel(testCase,level,side)
            import matlab.unittest.constraints.HasElementCount
            testCase.verifyThat(sierpinski(level), ...
                HasElementCount(side^2))
        end
    end 
end


