(omega-user-time-stepping)=

# Time stepping
Time stepper refers to the numerical scheme used to advance the model in time.
Omega implements a couple of time-stepping schemes. The user can select the
scheme they want in the configuration file.
```yaml
    TimeStepping:
       TimeStepperType: 'Forward-Backward'
```

The following time steppers are currently available:
| Config option name | Scheme |
| ------------------- | ------- |
| Forward-Backward | forward-backward |
| RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
