from pysb_Boolean import SimStep


def get_sim_steps(which):

    _supported = ["untreated", "30_min_BRAFi", "8_days_BRAFi", "8_days_BRAFi_plus_MEKi"]

    if which not in _supported:
        print("Error in 'pysb_Boolean::get_sim_steps': requested protocol ('%s') not found. Please try again." %
              str(which))
        quit()

    step_labels = None
    delta_ts = None
    conditions = None

    # Untreated
    if which == "untreated":

        step_labels = [
            "equilibration",
            "remove external calcium",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 50, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 30 min BRAFi
    elif which == "30_min_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 50, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi
    elif which == "8_days_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 500, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi + 15 min MEKi
    elif which == "8_days_BRAFi_plus_MEKi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add MEKi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 500, 10, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("MEKi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    return [SimStep(label, dt, cond) for label, dt, cond in zip(step_labels, delta_ts, conditions)]