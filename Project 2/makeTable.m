table1 = true;

% Part 1 and 2 table
if table1
    %Species for table



    %Table
    fprintf('\n')
    fprintf('Gas                      LHV(J/kg)       HHV(J/kg)        Exergy(J/kg)   Flow Exergy(J/kg)\n')
    fprintf('----------------------------------------------------------------------------------------------\n')
    fprintf('Hydrogen               %10.3e         %10.3e        %10.3e        %10.3e\n',LHV_values(1),HHV_values(1),int_ex_values(1),int_ex_flow_values(1))
    fprintf('Carbon Monoxide        %10.3e         %10.3e        %10.3e        %10.3e\n',LHV_values(2),HHV_values(2),int_ex_values(2),int_ex_flow_values(2))
    fprintf('Methane                %10.3e         %10.3e        %10.3e        %10.3e\n',LHV_values(3),HHV_values(3),int_ex_values(3),int_ex_flow_values(3))
    fprintf('Propane                %10.3e         %10.3e        %10.3e        %10.3e\n',LHV_values(4),HHV_values(4),int_ex_values(4),int_ex_flow_values(4))
    fprintf('Nitrogen               %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(5),int_ex_flow_values(5))
    fprintf('Oxygen                 %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(6),int_ex_flow_values(6))
    fprintf('Carbon Dioxide         %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,4.410e5,4.410e5)
    fprintf('Natural Gas            %10.3e         %10.3e        %10.3e        %10.3e\n',4.612e7,5.107e7,int_ex_values(8),int_ex_flow_values(8))
    fprintf('Simple Syngas          %10.3e         %10.3e        %10.3e        %10.3e\n',2.081e7,2.292e7,int_ex_values(9),int_ex_flow_values(9))
    fprintf('Engineering Air        %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(10),int_ex_flow_values(10))
    fprintf('Compressed Eng. Air    %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(11),int_ex_flow_values(11))
    fprintf('Cold Eng. Air          %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(12),int_ex_flow_values(12))
    fprintf('Warm Eng. Air          %10.3e         %10.3e        %10.3e        %10.3e\n',0,0,int_ex_values(13),int_ex_flow_values(13))
end