#' @title Generate experiment design
#'
#' @param n number of samples for the design
#' @param br number of biological samples for batch replicates
#' @param ncol number of columns per plate
#' @param nrow number of samples per plate
#'
#' @return A matrix of experiment design
#'
#' @export expDesign
expDesign = function(n = 1000, br = 5, ncol = 10, nrow = 10) {
    first_row = c("Pool 1", 1:(ncol-1))
    exp_mat = matrix(first_row, nrow = 1)

    cur_row_id = 2
    cur_sample = ncol
    tot_samples = ncol-1
    new_plate = FALSE

    num_p_plate = (ncol-1)*nrow
    # loop each row
    while(TRUE) {
        cur_row = NULL
        if (new_plate) { # Check if new plate
            result = get_new_plate_row(exp_mat, cur_row_id, cur_sample, tot_samples, ncol, br, plate_rows = nrow, new_plate, num_p_plate)
        } else {
            result = add_normal_row(exp_mat, cur_row, cur_row_id, tot_samples, ncol, cur_sample, new_plate, num_p_plate)
        }
        exp_mat = result$exp_mat
        tot_samples = result$tot_samples
        cur_sample = result$cur_sample
        new_plate = result$new_plate
        cur_row_id = increment_sample(cur_row_id)
        if (cur_sample > n) {
            last_row = exp_mat[nrow(exp_mat),]
            rmID = NULL
            if (any(grepl("\\*", last_row))) {
                rmID = which(grepl("\\*", last_row))
            }
            rmID = union(rmID, which(grepl("Pool", last_row)))
            # print(rmID)
            rm_samples = last_row[-rmID]
            rm_samples = as.numeric(rm_samples)
            if (any(rm_samples > n)) {
                rm_samples = rm_samples[which(rm_samples > n)]
                last_row = last_row[-which(last_row %in% rm_samples)]
            }
            last_row = c(last_row, rep(NA, ncol-length(last_row)))
            exp_mat[nrow(exp_mat),] = last_row
            break
        }
    }
    colnames(exp_mat) = c("Pooled sample", paste("col ", 1:(ncol-1), sep = ""))
    num_plates = ceiling(nrow(exp_mat)/nrow)
    num_sr = paste("SR", 1:nrow, sep = "")
    row_names = paste(rep(paste("Plate", 1:num_plates, sep = ""), each = length(num_sr)), "_", num_sr, sep = "")
    rownames(exp_mat) = row_names[1:nrow(exp_mat)]


    exp_mat

}




increment_sample = function(any_num) {
    any_num + 1
}

add_control = function(exp_mat, prev_row_id, long_control = FALSE) {
    prev_samples = as.character(exp_mat[prev_row_id,])
    prev_samples = prev_samples[-which(grepl("Pool", prev_samples))]
    if (any(grepl("\\*", prev_samples))) {
        prev_samples = prev_samples[-which(grepl("\\*", prev_samples))]
    }
    random_repeat = sample(prev_samples, 1)
    if (long_control) {
        random_repeat = paste(random_repeat, "**", sep = "")
    } else {
        random_repeat = paste(random_repeat, "*", sep = "")
    }
    random_repeat
}

add_long_control_row = function(exp_mat, cur_row, cur_row_id, tot_samples, num_p_sr, cur_sample, long_run, new_plate, num_p_plate) {
    if (is.null(cur_row)) { # not plate repeats
        next_id = 1
    } else {
        next_id = length(cur_row)+1
    }
    possible_id = num_p_sr:next_id
    repeat_id = sample(possible_id, 2)
    for (i in next_id:num_p_sr) {
        if (i == repeat_id[1]) { # long run control
            random_repeat = add_control(exp_mat, prev_row_id = cur_row_id-long_run+1, TRUE)
            cur_row = c(cur_row, random_repeat)
        } else if (i == repeat_id[2]) { # random control
            random_repeat = add_control(exp_mat, cur_row_id-1, FALSE)
            cur_row = c(cur_row, random_repeat)
        } else { # sample
            cur_row = c(cur_row, cur_sample)
            cur_sample = increment_sample(cur_sample)
        }
        tot_samples = increment_sample(tot_samples)
        if ((tot_samples %% num_p_plate) == 0) {
            new_plate = TRUE
            cur_row = c(cur_row, rep(NA, num_p_sr-length(cur_row)))
            break
        }
    }

    cur_row = c(paste("Pool ", cur_row_id, sep = ""), cur_row)
    exp_mat= rbind(exp_mat, cur_row)
    return(
        list(
            exp_mat = exp_mat,
            cur_sample = cur_sample,
            tot_samples = tot_samples,
            new_plate = new_plate
        )
    )
}

add_normal_row = function(exp_mat, cur_row, cur_row_id, tot_samples, num_p_sr, cur_sample, new_plate, num_p_plate) {
    if (is.null(cur_row)) {
        next_id = 1
    } else {
        next_id = length(cur_row)+1
    }
    num_p_sr = num_p_sr-1
    possible_id = num_p_sr:next_id
    repeat_id = sample(possible_id, 1)
    for (i in next_id:num_p_sr) {
        if (i == repeat_id) { # random control
            random_repeat = add_control(exp_mat, cur_row_id-1, FALSE)
            cur_row = c(cur_row, random_repeat)
        } else { # sample
            cur_row = c(cur_row, cur_sample)
            cur_sample = increment_sample(cur_sample)
        }
        tot_samples = increment_sample(tot_samples)
        if ((tot_samples %% num_p_plate) == 0) {
            new_plate = TRUE
            cur_row = c(cur_row, rep(NA, num_p_sr-length(cur_row)))
            break
        }
    }
    cur_row = c(paste("Pool ", cur_row_id, sep = ""), cur_row)
    exp_mat= rbind(exp_mat, cur_row)
    return(
        list(
            exp_mat = exp_mat,
            cur_sample = cur_sample,
            tot_samples = tot_samples,
            new_plate = new_plate
        )
    )
}

get_new_plate_row = function(exp_mat, cur_row_id, cur_sample, tot_samples, num_p_sr, plate_repeat_size, plate_rows, new_plate, num_p_plate) {
    prev_row_id = cur_row_id-1
    prev_plate = as.character(exp_mat[prev_row_id:(cur_row_id-plate_rows),])
    prev_plate = prev_plate[-which(grepl("Pool", prev_plate))]
    prev_plate = prev_plate[-which(grepl("\\*", prev_plate))]
    plate_repeats = sample(prev_plate, plate_repeat_size)
    plate_repeats = paste(plate_repeats, "***", sep = "")
    cur_row = plate_repeats
    for(i in 1:plate_repeat_size) {
        tot_samples = increment_sample(tot_samples)
    }
    # num_p_sr = num_p_sr-1

    if (length(cur_row) < (num_p_sr-1)) {
        # Check if long run repeat is necessary
        # if ((cur_row_id %% long_run) == 0) {
        #     # Add long control and go next
        #     result = add_long_control_row(exp_mat, cur_row, cur_row_id, tot_samples, num_p_sr, cur_sample, long_run, new_plate, num_p_plate)
        # } else {
        # Add regular samples (with repeats)
        result = add_normal_row(exp_mat, cur_row, cur_row_id, tot_samples, num_p_sr, cur_sample, new_plate, num_p_plate)
        # }
        exp_mat = result$exp_mat
        tot_samples = result$tot_samples
        cur_sample = result$cur_sample
        new_plate = FALSE
        # cur_row_id = increment_sample(cur_row_id)
    } else {
        stop("Too much plate repeats")
    }
    return(
        list(
            exp_mat = exp_mat,
            cur_sample = cur_sample,
            tot_samples = tot_samples,
            new_plate = new_plate
        )
    )
}
