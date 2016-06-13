function B = edge_mask(m,n,order,vertical,circ_lr,circ_tb)
% Computes an edge mask, for Jv (vertical) or Jh (horizontal) edges. Can be
% circular in left/right direction, but not in top/bottom direction (this
% reflects the current implementation of the algorithm).

if(~mod(order,2))
    error(message('smoothbox::evenOrder'));
end

if(~circ_lr && ~circ_tb) % Non-circular
    if(vertical)
        B = ones(m-1,n);
        K = ones(order)/(order^2);
        B = conv2(B,K,'same');
        B = [B ; zeros(1, n)];
    else
        B = ones(m,n-1);
        K = ones(order)/(order^2);
        B = conv2(B,K,'same');
        B = [B zeros(m,1)];
    end        
elseif(circ_lr && ~circ_tb) % Circular in left/right direction
    if(vertical)
        B = ones(m-1,1);
        K = ones(order,1)/order;
        B = conv(B,K,'same');
        B = [B ; 0];
        B = repmat(B, [1 n]);
    else
        B = ones(m,1);
        K = ones(order,1)/order;
        B = conv(B,K,'same');
        B = repmat(B, [1 n]);
    end
elseif(~circ_lr && circ_tb) % Circular in top/bottom direction
    if(vertical)
        B = ones(1,n);
        K = ones(1,order)/order;
        B = conv(B,K,'same');
        B = repmat(B, [m 1]);
    else
        B = ones(1,n-1);
        K = ones(1,order)/order;
        B = conv(B,K,'same');
        B = [B 0];
        B = repmat(B, [m 1]);
    end
else % Circular in both directions
    B = ones(m,n);
end