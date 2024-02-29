function print_monophyleticRobustnessProgress(iterspassed::Int64, totaliters::Int64, start_time::Float64)
    percent_complete = round(100 * iterspassed / totaliters, digits = 2)
    time_passed = time() - start_time
    est_total_runtime = time_passed / (percent_complete / 100)
    units = "seconds"

    # if more than 90 min passed, print in hours
    # elseif more than 180 sec passed, print in minutes
    if est_total_runtime > 60 * 90
        time_passed /= 60 * 60
        est_total_runtime /= 60 * 60
        units = "hours"
    elseif est_total_runtime > 180
        time_passed /= 60
        est_total_runtime /= 60
        units = "minutes"
    end
    time_passed = round(time_passed, digits = 2)
    est_total_runtime = round(est_total_runtime, digits = 2)

    print("\r\t$(percent_complete)% ($(iterspassed)/$(totaliters)) complete "*
          "($(time_passed)/≈$(est_total_runtime) $(units))              ")
end