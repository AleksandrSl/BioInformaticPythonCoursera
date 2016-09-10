# Calculates the minimum number of coins with denominations from "coins" in the change

# !!!!!Recall that our original goal was to make change, not just compute MinNumCoins(money). Modify DPChange so that
# it not only computes the minimum number of coins but also returns these coins.

def dynamic_progr_change(money, coins):

    min_num_of_coins = [0]
    for m in range(1, money + 1):
        if m > max(coins):
            min_num_of_coins.pop(0)
        min_num_of_coins.append(10000)
        for coin in coins:
            if m > max(coins):
                previous_min = min_num_of_coins[-coin + 1] + 1
                if previous_min < min_num_of_coins[-1]:
                    min_num_of_coins[-1] = previous_min
            else:
                if coin > m:
                    continue
                previous_min = min_num_of_coins[(m - coin)] + 1
                if previous_min < min_num_of_coins[m]:
                    min_num_of_coins[m] = previous_min
    return min_num_of_coins[-1]

print(dynamic_progr_change(19418, [23, 5, 3, 1]))