function plot_trajectories(t,x)

% ---------- FORCE LIGHT MODE (GLOBAL) ----------
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');

% Time is always 0 at the beginning of an ODE45 solution
ind0 = find(t == 0);
m = length(ind0);

ind0(length(ind0)+1) = length(t)+1;

t_tot = [];

for i = 1:m
    
    if isempty(t_tot)
        t_end = 0;
    else
        t_end = t_tot(end);
    end
    
    t_tot = [t_tot; t_end + t(ind0(i):ind0(i+1)-1)];
    
end

%% ------------------ Figure 1: Phase Portrait ------------------ %%
figure(1)
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w', 'XColor','k','YColor','k');

plot(x(:,1), x(:,4))
grid on
title('Phase portrait of q_1 vs dq_1')
xlabel('q_1 (rad)')
ylabel('dq_1 (rad/s)')

exportgraphics(gcf, 'phase_portrait_full_dynamics.pdf', ...
    'ContentType','vector', ...
    'Resolution',300);


%% ------------------ Figure 2: Joint Angles ------------------ %%
figure(2)
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w', 'XColor','k','YColor','k');

plot(t_tot, x(:,1)), hold on
plot(t_tot, x(:,2))
plot(t_tot, x(:,3))
hold off, grid on

title('Joint angles')
xlabel('t (s)')
ylabel('q (rad)')
legend('q_1','q_2','q_3')

exportgraphics(gcf, 'joint_angles_full_dynamics.pdf', ...
    'ContentType','vector', ...
    'Resolution',300);


%% ------------------ Figure 3: Joint Velocities ------------------ %%
figure(3)
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w', 'XColor','k','YColor','k');

plot(t_tot, x(:,4)), hold on
plot(t_tot, x(:,5))
plot(t_tot, x(:,6))
hold off, grid on

title('Joint velocities')
xlabel('t (s)')
ylabel('dq (rad/s)')
legend('dq_1','dq_2','dq_3')

exportgraphics(gcf, 'joint_velocities_full_dynamics.pdf', ...
    'ContentType','vector', ...
    'Resolution',300);

end