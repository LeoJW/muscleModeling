% MINIMIZATION ON COST FOR TPB VS EMR MANUS EXTENSION
clear variables; 


%% Define Constants
niter = 10;
niiter = 10;

k = -1/3;
mbody = linspace(0.005,1,niter); %[kg]

ltpb = 0.05; %0-1
lemr = 0.25; %0-1
    
tpbAngle = linspace(5,25,niiter); %degrees

%% Iterate and Minimize

%---Define Constraints (Note tpb will always be 1st entry)
%-Equality constraints
Aeq = [sind(tpbAngle),1];
beq = mbody*9.81*0.25; %(Required manus force)
% beq = ones(size(mbody));
%-Bounds
lb = [zeros(niter,1),0.25*beq.'];
ub = [10,10];

%Preallocate
M = zeros(2,niter,niiter);
x0 = zeros(2,niter,niiter);

for ii = 1:niiter
    
    Aeq = [sind(tpbAngle(ii)),1];

    %Loop thru body masses
    for i = 1:niter
        %Define objective function
        fcost = @(masses) mbody(i)^(2*k)*(masses(1)*ltpb + masses(2)*lemr);
        %Define intial conditions
        x0(2,i,ii) = rand*beq(i); %EMR
        x0(1,i,ii) = (beq(i) - x0(2,i,ii))/sind(tpbAngle(ii)); %TPB

        %Minimize!
        [M(:,i,ii)] = fmincon(fcost,x0(:,i,ii),[],[],Aeq,beq(i),lb(i,:),ub);
    end
end


%% Plot Results


tpbcols = spring(niiter);
emrcols = winter(niiter);


figure()
hold on
grid on

for ii = 1:niiter
%     plot(mbody,M(1,:,ii)*sind(tpbAngle(ii))/beq,'.','color',emrcols(ii,:))
%     plot(mbody,M(2,:,ii)/beq,'o','color',emrcols(ii,:))
    plot(mbody,M(1,:,ii)*sind(tpbAngle(ii))/beq,'.','color',emrcols(ii,:))
    plot(mbody,M(2,:,ii)/beq,'o','color',emrcols(ii,:))
end

xlabel('Mbody [kg]')
% ylabel('Prop. of Manus Force')


tpbAngle(M(1,1,:)>1e-3)

%% Notes

%---Seems always to go "all-or-none": Either the TPB gets all, or the EMR
%gets all. Will flip between the two when tpbAngle is sufficiently high.
%For ltpb=0.09 and lemr=0.25, this is around 20-21 degrees
%---The tpbAngle at which optimal muscle flips is likely set by the balance
%between ltpb & lemr