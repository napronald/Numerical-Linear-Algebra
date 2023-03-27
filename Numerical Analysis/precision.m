% 6-bit system
clc
% for n1 = 0:1
%     for n2 = 0:1
%         for n3 = 0:1
%             for n4 = 0:1
%                 for n5 = 0:1
%                     for n6 = 0:1
%                         input = [n1 n2 n3 n4 n5 n6]
%                         if input(1) == 0
%                             sign = 1;
%                         else 
%                             sign = -1;
%                         end
%                         c = input(2)*2^1+input(3)*2^0;
%                         f = input(4)*(1/2)+input(5)*(1/4)+input(6)*(1/8);
%                         number_6bit = sign*2^(c-1)*(1+f)
%                         figure(1)
%                         hold on
%                         plot(number_6bit, 0,'bo')
%                     end
%                 end
%             end
%         end
%     end
% end

% %8-bit system
for n1 = 0:1
    for n2 = 0:1
        for n3 = 0:1
            for n4 = 0:1
                for n5 = 0:1
                    for n6 = 0:1
                        for n7 = 0:1
                            for n8 = 0:1
                                for n9 = 0:1
                                    input = [n1 n2 n3 n4 n5 n6 n7 n8 n9]
                                    if input(1) == 0
                                        sign = 1;
                                    else 
                                        sign = -1;
                                    end
                                    c = input(2)*2^3 + input(3)*2^2 + input(4)*2^1 + input(5)*2^0;
                                    f = input(6)*(1/2) + input(7)*(1/4) + input(8)*(1/8) + input(9)*(1/16);
                                    number_8bit = sign*2^(c-7)*(1+f)
                                    figure(2)
                                    hold on
                                    plot(number_8bit, 0,'bo')
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% What can you say about distribution of these numbers on the real line in 
% terms of how well they cover the represented interval?
% ANSWER:The numbers represented will be determined by the normalizer set.
% The normalizer set in my code is 2 making the lowest nonzero value 0.25
% and 3.75 the largest number. Numbers near 0.25 are more precise compared 
% to values near 3.75

% What is the largest and what is the smallest in absolute value number 
% that can be represented with this system?
% ANSWER:With a normalizer of 2, the largest number that can be represented 
% is 3.750, the smallest absolute value number that can be represented is 
% 0.250

% If all numbers on the interval between the smallest and the largest 
% numbers are represented with this 6-bit number system, which parts of 
% this interval will numbers represented with the smallest absolute error 
% in the representation and which ones have smallest relative error? 
% Give examples.
% ANSWER:Representing smaller numbers will have smaller absolute and
% relative error whereas larger numbers will have larger absolute and
% relative errors. Example 1. p=0.4 p*=0.40625 abs error = 0.0625 rel 
% error 0.015625 Example 2 p=3 p*=2.85 abs error=0.15 rel error = 0.06666

% Modify your program to increase the bits available to 8 and plot the 
% numbers represented with 4 bits for the exponent and 4 bits for the 
% mantissa. What changed?
% ANSWER:We have more numbers that can be represented making the accurracy 
% of our approximation stronger for small numbers but worse for big numbers

% Are the small numbers between 0 and 1 represented well by the 6- and 
% 8-bit systems? If yes, which of these systems represents the interval 
% [0,1] better? If no, then how can you modify your finite-precision system
% to get this interval to be represented better?
% ANSWER:Of the two systems, the 8-bit system better represents small 
% numbers on the interval [0,1]. Currently, small numbers between 0 and 1 
% could be better represented. We can improve both systems by adjusting our
% normalizer to a larger number.



