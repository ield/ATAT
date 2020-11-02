function [sll] = calcSLL(y, plane)
    % Function created tp find the sll. It is found the second maximum of
    % the function and it is returned its value
    limDer = 0.1;   % Limit impossed to determine the places where there is maximum or minimum. (since d is not 0 never because of discretization).
    ydiff = diff(y(10:end));
    maxim = ydiff < limDer & ydiff > -limDer;
    index = find(maxim == 1, 2);
    sll = y(index(2)+10 + 1);
    fprintf('The SLL in %s are %f dB\n', plane, sll);
end

