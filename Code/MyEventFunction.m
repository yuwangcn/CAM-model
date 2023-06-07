function [VALUE, ISTERMINAL, DIRECTION] = MyEventFunction(~,~)

    TimeOut = 60; % 20min
%     toc
    VALUE = toc-TimeOut;
    ISTERMINAL=1;
    DIRECTION = 0;
    if VALUE>0
        errorStruct.message='Time out. Possibily singular values';
        errorStruct.identifier='eCM:Timeout';
        error(errorStruct);
    end
end