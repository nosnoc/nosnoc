%% experimetns


% with mpcc mode 5 as usual
try
    main_oscilator_integration_order_via_M
catch
    disp('oh')
end


try
    main_oscilator_integration_order_via_M_fesd
catch
    disp('oh')
end

%% with mppc mode 2 (sanity check)
try
    main_oscilator_integration_order
catch
    disp('oh')
end


try
    main_oscilator_integration_order_fesd
catch
    disp('oh')
end